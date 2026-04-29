%% Tutorial 06: Periodic Ewald solver interface
%
% This tutorial demonstrates the periodic-Ewald analogue of Tutorial 05.
%
% Goal:
%   Build a small periodic polarizable system, prepare an SCF problem, build
%   solver-specific periodic Ewald polarization operators, solve with the
%   available solvers, and compare the results.
%
% This tutorial exercises:
%   - thole.prepare_scf_problem
%   - thole.make_polarization_operator(..., 'Mode','periodic_ewald')
%   - thole.solve_scf_direct
%   - thole.solve_scf_jacobi
%   - thole.solve_scf_gmres
%   - thole.solve_scf_sor
%   - calc.compute_total_energy_active_space
%
% Periodic conventions:
%   Polarize uses direct lattice vectors as ROWS:
%
%       cart = frac * H
%
%   The periodic Ewald routines use geom.get_lattice / 
%   ewald.enumerate_kvecs_from_lattice, so public periodic code should not
%   transpose the lattice manually.

clear;
clc;

fprintf('\n============================================================\n');
fprintf('Tutorial 06: Periodic Ewald solver interface\n');
fprintf('============================================================\n\n');

%% ------------------------------------------------------------------------
% 1. Build a tiny periodic polarizable system
% -------------------------------------------------------------------------

sys = local_make_periodic_toy_system();

fprintf('[1] Periodic toy system\n');
fprintf('  n_sites        = %d\n', sys.n_sites);
fprintf('  n polarizable  = %d\n', nnz(sys.site_is_polarizable));
fprintf('  is_periodic    = %d\n', sys.is_periodic);
fprintf('  lattice rows H =\n');
disp(sys.super_lattice);

lat = geom.get_lattice(sys);

fprintf('  lattice volume = %.8f bohr^3\n', lat.volume);
fprintf('  ||H*G - 2*pi*I||_F = %.3e\n\n', ...
    norm(lat.H * lat.G - 2*pi*eye(3), 'fro'));

io.assert_atomic_units(sys);

%% ------------------------------------------------------------------------
% 2. Define an external field on active polarizable sites
% -------------------------------------------------------------------------
%
% In a real workflow, this Eext would often come from:
%
%   calc.compute_external_field(...)
%
% for a charged molecule/pair embedded in a polarizable periodic crystal.
%
% Here we set a small explicit field so the tutorial focuses on the solver
% interface and periodic induced-dipole operator.

Eext = zeros(sys.n_sites, 3);

Eext(1, :) = [ 1.0e-3, -0.5e-3,  0.2e-3];
Eext(2, :) = [-0.4e-3,  0.8e-3, -0.1e-3];
Eext(3, :) = [ 0.3e-3,  0.1e-3,  0.5e-3];

fprintf('[2] External field\n');
fprintf('  ||Eext||_F = %.8e\n\n', norm(Eext, 'fro'));

%% ------------------------------------------------------------------------
% 3. Prepare SCF problem
% -------------------------------------------------------------------------

scfParams = struct();
scfParams.tol = 1e-11;
scfParams.maxIter = 500;
scfParams.mixing = 0.4;
scfParams.omega = 1.0;
scfParams.verbose = false;

problem = thole.prepare_scf_problem(sys, Eext, scfParams);

fprintf('[3] SCF problem\n');
fprintf('  nSites        = %d\n', problem.nSites);
fprintf('  nPolSites     = %d\n', problem.nPolSites);
fprintf('  activeSites   = [%s]\n', sprintf('%d ', problem.activeSites));
fprintf('  ||Eext_pol||  = %.8e\n', norm(problem.Eext_pol_vec));
fprintf('  ||alpha_vec|| = %.8e\n\n', norm(problem.alpha_pol_vec));

%% ------------------------------------------------------------------------
% 4. Periodic Ewald operator parameters
% -------------------------------------------------------------------------
%
% These parameters are deliberately small because this is a tiny tutorial.
% For production calculations, rcut/kcut/alpha should be selected and
% converged carefully.
%
% Important:
%   geom.build_active_row_cache_periodic enforces rcut < Lmin/2.
%
% In this toy lattice, Lmin is comfortably larger than the chosen rcut.

ewaldParams = struct();
ewaldParams.alpha = 0.30;
ewaldParams.rcut = 9.0;
ewaldParams.kcut = 1.25;
ewaldParams.boundary = 'tinfoil';

Lmin = geom.shortest_lattice_translation(sys.super_lattice);

fprintf('[4] Periodic Ewald parameters\n');
fprintf('  alpha     = %.6f\n', ewaldParams.alpha);
fprintf('  rcut      = %.6f bohr\n', ewaldParams.rcut);
fprintf('  kcut      = %.6f bohr^-1\n', ewaldParams.kcut);
fprintf('  boundary  = %s\n', ewaldParams.boundary);
fprintf('  Lmin      = %.6f bohr\n', Lmin);
fprintf('  Lmin/2    = %.6f bohr\n', 0.5 * Lmin);

if ~(ewaldParams.rcut < 0.5 * Lmin)
    error('Tutorial rcut violates the single-image condition rcut < Lmin/2.');
end

fprintf('\n');

baseOpArgs = {
    'Mode', 'periodic_ewald', ...
    'UseThole', true, ...
    'Rcut', ewaldParams.rcut, ...
    'Alpha', ewaldParams.alpha, ...
    'Kcut', ewaldParams.kcut, ...
    'Boundary', ewaldParams.boundary, ...
    'KspaceMode', 'full', ...
    'UseMex', true, ...
    'Profile', false, ...
    'Verbose', false ...
};

%% ------------------------------------------------------------------------
% 5. Direct solver / dense periodic operator
% -------------------------------------------------------------------------
%
% The dense periodic operator is intended for tiny validation systems only.
% It is useful here as a reference solution.

fprintf('[5] Direct solve with dense periodic operator\n');

tBuild = tic;
opDirect = thole.make_polarization_operator( ...
    sys, problem, baseOpArgs{:}, ...
    'Solver', 'direct', ...
    'Backend', 'auto');
buildDirectTime = toc(tBuild);

fprintf('  op.kind       = %s\n', opDirect.kind);
fprintf('  op.backend    = %s\n', opDirect.backend);
fprintf('  op.size       = [%d %d]\n', opDirect.size(1), opDirect.size(2));
fprintf('  build time    = %.6f s\n', buildDirectTime);
fprintf('  nK            = %d\n', opDirect.info.nK);
fprintf('  nRealEntries  = %d\n', opDirect.info.nRealEntriesDirected);

tSolve = tic;
[muDirect, infoDirect] = thole.solve_scf_direct(problem, opDirect);
solveDirectTime = toc(tSolve);

fprintf('  converged     = %d\n', infoDirect.converged);
fprintf('  relres        = %.3e\n', infoDirect.relres);
fprintf('  solve time    = %.6f s\n', solveDirectTime);
fprintf('  ||mu||_F      = %.8e\n\n', norm(muDirect, 'fro'));

%% ------------------------------------------------------------------------
% 6. Jacobi solver / periodic paircache operator
% -------------------------------------------------------------------------

fprintf('[6] Jacobi solve with periodic paircache operator\n');

tBuild = tic;
opJacobi = thole.make_polarization_operator( ...
    sys, problem, baseOpArgs{:}, ...
    'Solver', 'jacobi', ...
    'Backend', 'auto');
buildJacobiTime = toc(tBuild);

fprintf('  op.kind       = %s\n', opJacobi.kind);
fprintf('  op.backend    = %s\n', opJacobi.backend);
fprintf('  row_update    = %d\n', opJacobi.capabilities.row_update);
fprintf('  build time    = %.6f s\n', buildJacobiTime);
fprintf('  nK            = %d\n', opJacobi.info.nK);
fprintf('  nRealEntries  = %d\n', opJacobi.info.nRealEntriesDirected);

jacobiOpts = struct();
jacobiOpts.tol = 1e-10;
jacobiOpts.max_iter = 1000;
jacobiOpts.mixing = 0.35;
jacobiOpts.stop_metric = 'relres';
jacobiOpts.verbose = false;

tSolve = tic;
[muJacobi, infoJacobi] = thole.solve_scf_jacobi(problem, opJacobi, jacobiOpts);
solveJacobiTime = toc(tSolve);

fprintf('  converged     = %d\n', infoJacobi.converged);
fprintf('  iterations    = %d\n', infoJacobi.iterations);
fprintf('  relres        = %.3e\n', infoJacobi.relres);
fprintf('  solve time    = %.6f s\n', solveJacobiTime);
fprintf('  ||mu||_F      = %.8e\n', norm(muJacobi, 'fro'));
fprintf('  vs direct     = %.3e\n\n', local_active_mu_error(muJacobi, muDirect, problem));

%% ------------------------------------------------------------------------
% 7. GMRES solver / periodic paircache operator
% -------------------------------------------------------------------------

fprintf('[7] GMRES solve with periodic paircache operator\n');

tBuild = tic;
opGmres = thole.make_polarization_operator( ...
    sys, problem, baseOpArgs{:}, ...
    'Solver', 'gmres', ...
    'Backend', 'auto');
buildGmresTime = toc(tBuild);

fprintf('  op.kind       = %s\n', opGmres.kind);
fprintf('  op.backend    = %s\n', opGmres.backend);
fprintf('  row_update    = %d\n', opGmres.capabilities.row_update);
fprintf('  build time    = %.6f s\n', buildGmresTime);
fprintf('  nK            = %d\n', opGmres.info.nK);
fprintf('  nRealEntries  = %d\n', opGmres.info.nRealEntriesDirected);

gmresOpts = struct();
gmresOpts.tol = 1e-11;
gmresOpts.max_iter = 100;
gmresOpts.verbose = false;

tSolve = tic;
[muGmres, infoGmres] = thole.solve_scf_gmres(problem, opGmres, gmresOpts);
solveGmresTime = toc(tSolve);

fprintf('  converged     = %d\n', infoGmres.converged);
fprintf('  iterations    = %d\n', infoGmres.iterations);
fprintf('  relres        = %.3e\n', infoGmres.relres);
fprintf('  solve time    = %.6f s\n', solveGmresTime);
fprintf('  ||mu||_F      = %.8e\n', norm(muGmres, 'fro'));
fprintf('  vs direct     = %.3e\n\n', local_active_mu_error(muGmres, muDirect, problem));

%% ------------------------------------------------------------------------
% 8. SOR solver / periodic rowcache operator
% -------------------------------------------------------------------------
%
% This should use the periodic raw fast path inside thole.solve_scf_sor.
% That path directly contracts the real-space row cache and updates
% reciprocal A/B source sums incrementally during the sweep.

fprintf('[8] SOR solve with periodic rowcache operator\n');

tBuild = tic;
opSor = thole.make_polarization_operator( ...
    sys, problem, baseOpArgs{:}, ...
    'Solver', 'sor', ...
    'Backend', 'auto');
buildSorTime = toc(tBuild);

fprintf('  op.kind       = %s\n', opSor.kind);
fprintf('  op.backend    = %s\n', opSor.backend);
fprintf('  row_update    = %d\n', opSor.capabilities.row_update);
fprintf('  build time    = %.6f s\n', buildSorTime);
fprintf('  nK            = %d\n', opSor.info.nK);
fprintf('  nRealEntries  = %d\n', opSor.info.nRealEntriesDirected);

sorOpts = struct();
sorOpts.tol = 1e-10;
sorOpts.max_iter = 1000;
sorOpts.omega = 1.0;
sorOpts.stop_metric = 'relres';
sorOpts.residual_every = 1;
sorOpts.verbose = false;

tSolve = tic;
[muSor, infoSor] = thole.solve_scf_sor(problem, opSor, sorOpts);
solveSorTime = toc(tSolve);

fprintf('  converged     = %d\n', infoSor.converged);
fprintf('  iterations    = %d\n', infoSor.iterations);
fprintf('  relres        = %.3e\n', infoSor.relres);
fprintf('  solve time    = %.6f s\n', solveSorTime);
fprintf('  fast path     = %s\n', infoSor.rowcache_fast_path_type);
fprintf('  ||mu||_F      = %.8e\n', norm(muSor, 'fro'));
fprintf('  vs direct     = %.3e\n\n', local_active_mu_error(muSor, muDirect, problem));

%% ------------------------------------------------------------------------
% 9. Blocked k-space apply sanity check
% -------------------------------------------------------------------------
%
% This exercises the blocked reciprocal-space storage/apply path and checks
% it against full storage for one arbitrary active-space dipole vector.

fprintf('[9] Full vs blocked k-space apply check\n');

opBlocked = thole.make_polarization_operator( ...
    sys, problem, baseOpArgs{:}, ...
    'Solver', 'gmres', ...
    'Backend', 'auto', ...
    'KspaceMode', 'blocked', ...
    'KBlockSize', 3);

testMu = zeros(3 * problem.nPolSites, 1);
testMu(1:3:end) = [ 1.0e-4; -2.0e-4;  0.5e-4];
testMu(2:3:end) = [ 2.0e-4;  1.5e-4; -1.0e-4];
testMu(3:3:end) = [-1.0e-4;  0.5e-4;  2.0e-4];

yFull = opGmres.apply(testMu);
yBlocked = opBlocked.apply(testMu);

blockedErr = norm(yFull - yBlocked) / max(1, norm(yFull));

fprintf('  full storage mode     = %s\n', opGmres.k_cache.storage_mode);
fprintf('  blocked storage mode  = %s\n', opBlocked.k_cache.storage_mode);
fprintf('  blocked nBlocks       = %d\n', opBlocked.k_cache.num_blocks);
fprintf('  relative apply error  = %.3e\n\n', blockedErr);

if blockedErr > 1e-12
    error('Blocked periodic apply path does not match full path.');
end

%% ------------------------------------------------------------------------
% 10. Energy comparison
% -------------------------------------------------------------------------
%
% calc.compute_total_energy_active_space uses the stationary active-space
% SCF expression and can use op.apply for matrix-free operators.

fprintf('[10] Active-space energy comparison\n');

Edir = calc.compute_total_energy_active_space(sys, problem, muDirect, Eext, opDirect);
Ejac = calc.compute_total_energy_active_space(sys, problem, muJacobi, Eext, opJacobi);
Egm  = calc.compute_total_energy_active_space(sys, problem, muGmres,  Eext, opGmres);
Esor = calc.compute_total_energy_active_space(sys, problem, muSor,    Eext, opSor);

fprintf('  direct total = %.12e Ha | stationary diff = %.3e Ha | relres = %.3e\n', ...
    Edir.total, Edir.stationary_consistency, Edir.relres);

fprintf('  jacobi total = %.12e Ha | diff = %.3e Ha | stationary diff = %.3e Ha | relres = %.3e\n', ...
    Ejac.total, Ejac.total - Edir.total, Ejac.stationary_consistency, Ejac.relres);

fprintf('  gmres  total = %.12e Ha | diff = %.3e Ha | stationary diff = %.3e Ha | relres = %.3e\n', ...
    Egm.total, Egm.total - Edir.total, Egm.stationary_consistency, Egm.relres);

fprintf('  sor    total = %.12e Ha | diff = %.3e Ha | stationary diff = %.3e Ha | relres = %.3e\n\n', ...
    Esor.total, Esor.total - Edir.total, Esor.stationary_consistency, Esor.relres);

%% ------------------------------------------------------------------------
% 11. Summary
% -------------------------------------------------------------------------

fprintf('============================================================\n');
fprintf('Periodic Ewald solver summary\n');
fprintf('============================================================\n');

fprintf('%-10s %-26s %-10s %-12s %-12s %-12s\n', ...
    'solver', 'backend', 'converged', 'iters', 'relres', 'vs direct');

fprintf('%-10s %-26s %-10d %-12s %-12.3e %-12.3e\n', ...
    'direct', opDirect.backend, infoDirect.converged, '-', ...
    infoDirect.relres, 0.0);

fprintf('%-10s %-26s %-10d %-12d %-12.3e %-12.3e\n', ...
    'jacobi', opJacobi.backend, infoJacobi.converged, ...
    infoJacobi.iterations, infoJacobi.relres, ...
    local_active_mu_error(muJacobi, muDirect, problem));

fprintf('%-10s %-26s %-10d %-12d %-12.3e %-12.3e\n', ...
    'gmres', opGmres.backend, infoGmres.converged, ...
    infoGmres.iterations, infoGmres.relres, ...
    local_active_mu_error(muGmres, muDirect, problem));

fprintf('%-10s %-26s %-10d %-12d %-12.3e %-12.3e\n', ...
    'sor', opSor.backend, infoSor.converged, ...
    infoSor.iterations, infoSor.relres, ...
    local_active_mu_error(muSor, muDirect, problem));

fprintf('\n');
fprintf('Periodic SOR fast path: %s\n', infoSor.rowcache_fast_path_type);
fprintf('Done.\n');

%% =========================================================================
% Local helpers
% =========================================================================

function sys = local_make_periodic_toy_system()
sys = struct();

% Three active polarizable sites plus one nonpolarizable spectator.
%
% The mildly triclinic row-lattice prevents hidden row/column convention
% errors. Coordinates are Cartesian bohr.
sys.site_pos = [
     3.0   2.0   2.0
    10.0   6.0   5.0
    17.0  11.0   9.0
     6.0  15.0  12.0
];

sys.site_charge = zeros(4, 1);

% Small polarizabilities keep the tiny fixed-point problem comfortably
% stable for all solvers.
sys.site_alpha = [
    0.08
    0.07
    0.06
    0.00
];

sys.site_is_polarizable = [
    true
    true
    true
    false
];

sys.n_sites = 4;

sys.thole_a = 0.39;

sys.site_type = {'X'; 'X'; 'X'; 'X'};
sys.site_class = {'pol'; 'pol'; 'pol'; 'spectator'};
sys.site_label = {'p1'; 'p2'; 'p3'; 'n'};
sys.site_mol_id = [1; 2; 3; 4];
sys.site_is_active = sys.site_is_polarizable;

sys.units = struct();
sys.units.length = 'bohr';
sys.units.alpha = 'atomic_unit';
sys.units.charge = 'elementary_charge';

sys.is_periodic = true;
sys.periodic_mode = 'periodic';

% Row-lattice convention: cart = frac * H.
sys.lattice = [
    28.0   0.0   0.0
     2.0  27.0   0.0
     1.0   3.0  26.0
];

sys.super_lattice = sys.lattice;
end

function err = local_active_mu_error(mu, muRef, problem)
muAct = mu(problem.activeSites, :);
muRefAct = muRef(problem.activeSites, :);

err = norm(muAct - muRefAct, 'fro') / max(1, norm(muRefAct, 'fro'));
end