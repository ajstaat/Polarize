function test_thole_periodic_ewald_solvers()
%TEST_THOLE_PERIODIC_EWALD_SOLVERS Smoke/regression test for periodic Ewald solvers.
%
% Checks:
%   - periodic dense/direct operator builds
%   - periodic paircache operator builds for Jacobi and GMRES
%   - periodic rowcache operator builds for SOR
%   - all solver results agree with dense direct reference
%   - SOR uses the periodic raw row-cache fast path

sys = local_make_periodic_polsys();

Eext = zeros(sys.n_sites, 3);
Eext(1, :) = [ 1.0e-3, -0.5e-3,  0.2e-3];
Eext(2, :) = [-0.4e-3,  0.8e-3, -0.1e-3];
Eext(3, :) = [ 0.3e-3,  0.1e-3,  0.5e-3];

scfParams = struct();
scfParams.tol = 1e-11;
scfParams.maxIter = 500;
scfParams.mixing = 0.4;
scfParams.omega = 1.0;
scfParams.verbose = false;

problem = thole.prepare_scf_problem(sys, Eext, scfParams);

opArgs = {
    'Mode', 'periodic_ewald', ...
    'UseThole', true, ...
    'Rcut', 9.0, ...
    'Alpha', 0.30, ...
    'Kcut', 1.25, ...
    'Boundary', 'tinfoil', ...
    'KspaceMode', 'full', ...
    'UseMex', true, ...
    'Profile', false, ...
    'Verbose', false ...
};

% -------------------------------------------------------------------------
% Direct / dense reference.
% -------------------------------------------------------------------------
opDirect = thole.make_polarization_operator( ...
    sys, problem, opArgs{:}, ...
    'Solver', 'direct', ...
    'Backend', 'auto');

assert(strcmp(opDirect.mode, 'periodic_ewald'));
assert(strcmp(opDirect.kind, 'dense_matrix'));
assert(opDirect.capabilities.dense_matrix);

[muDirect, infoDirect] = thole.solve_scf_direct(problem, opDirect);

assert(all(isfinite(muDirect(:))), ...
    'Periodic direct solver returned non-finite dipoles.');

assert(infoDirect.relres < 1e-9, ...
    'Periodic direct solver residual is too large.');

muDirectPol = muDirect(problem.activeSites, :);

% -------------------------------------------------------------------------
% Jacobi / paircache.
% -------------------------------------------------------------------------
opJacobi = thole.make_polarization_operator( ...
    sys, problem, opArgs{:}, ...
    'Solver', 'jacobi', ...
    'Backend', 'auto');

assert(strcmp(opJacobi.mode, 'periodic_ewald'));
assert(strcmp(opJacobi.kind, 'matrix_free'));
assert(strcmp(opJacobi.backend, 'periodic_paircache_apply'));
assert(opJacobi.capabilities.apply);
assert(~opJacobi.capabilities.row_update);

jacobiOpts = struct();
jacobiOpts.tol = 1e-10;
jacobiOpts.max_iter = 1000;
jacobiOpts.mixing = 0.35;
jacobiOpts.stop_metric = 'relres';
jacobiOpts.verbose = false;

[muJacobi, infoJacobi] = thole.solve_scf_jacobi(problem, opJacobi, jacobiOpts);

assert(infoJacobi.converged, ...
    'Periodic Jacobi solver did not converge.');

assert(infoJacobi.relres < 1e-8, ...
    'Periodic Jacobi residual is too large.');

local_assert_mu_close(muJacobi(problem.activeSites, :), muDirectPol, ...
    5e-8, 'Jacobi periodic result differs from direct reference.');

% -------------------------------------------------------------------------
% GMRES / paircache.
% -------------------------------------------------------------------------
opGmres = thole.make_polarization_operator( ...
    sys, problem, opArgs{:}, ...
    'Solver', 'gmres', ...
    'Backend', 'auto');

assert(strcmp(opGmres.mode, 'periodic_ewald'));
assert(strcmp(opGmres.kind, 'matrix_free'));
assert(strcmp(opGmres.backend, 'periodic_paircache_apply'));
assert(opGmres.capabilities.apply);
assert(~opGmres.capabilities.row_update);

gmresOpts = struct();
gmresOpts.tol = 1e-11;
gmresOpts.max_iter = 100;
gmresOpts.verbose = false;

[muGmres, infoGmres] = thole.solve_scf_gmres(problem, opGmres, gmresOpts);

assert(infoGmres.converged, ...
    'Periodic GMRES solver did not converge.');

assert(infoGmres.relres < 1e-8, ...
    'Periodic GMRES residual is too large.');

local_assert_mu_close(muGmres(problem.activeSites, :), muDirectPol, ...
    5e-9, 'GMRES periodic result differs from direct reference.');

% -------------------------------------------------------------------------
% SOR / rowcache.
% -------------------------------------------------------------------------
opSor = thole.make_polarization_operator( ...
    sys, problem, opArgs{:}, ...
    'Solver', 'sor', ...
    'Backend', 'auto');

assert(strcmp(opSor.mode, 'periodic_ewald'));
assert(strcmp(opSor.kind, 'matrix_free'));
assert(strcmp(opSor.backend, 'periodic_rowcache_apply'));
assert(opSor.capabilities.apply);
assert(opSor.capabilities.row_update);

sorOpts = struct();
sorOpts.tol = 1e-10;
sorOpts.max_iter = 1000;
sorOpts.omega = 1.0;
sorOpts.stop_metric = 'relres';
sorOpts.residual_every = 1;
sorOpts.verbose = false;

[muSor, infoSor] = thole.solve_scf_sor(problem, opSor, sorOpts);

assert(infoSor.converged, ...
    'Periodic SOR solver did not converge.');

assert(infoSor.relres < 1e-8, ...
    'Periodic SOR residual is too large.');

assert(isfield(infoSor, 'used_periodic_fast_path') && infoSor.used_periodic_fast_path, ...
    'Periodic SOR should use the periodic raw fast path.');

assert(strcmp(infoSor.rowcache_fast_path_type, 'periodic_raw_rowcache'), ...
    'Periodic SOR should report periodic_raw_rowcache fast path.');

local_assert_mu_close(muSor(problem.activeSites, :), muDirectPol, ...
    5e-8, 'SOR periodic result differs from direct reference.');

% -------------------------------------------------------------------------
% Blocked k-space smoke test for the matrix-free apply path.
% -------------------------------------------------------------------------
opBlocked = thole.make_polarization_operator( ...
    sys, problem, opArgs{:}, ...
    'Solver', 'gmres', ...
    'Backend', 'auto', ...
    'KspaceMode', 'blocked', ...
    'KBlockSize', 3);

assert(strcmp(opBlocked.k_cache.storage_mode, 'blocked'), ...
    'Blocked periodic k-space operator should use blocked storage mode.');

testVec = util.stack_xyz([
     1.0e-4,  2.0e-4, -1.0e-4
    -2.0e-4,  1.5e-4,  0.5e-4
     0.5e-4, -1.0e-4,  2.0e-4
]);

yFull = opGmres.apply(testVec);
yBlocked = opBlocked.apply(testVec);

assert(norm(yFull - yBlocked) < 1e-12 * max(1, norm(yFull)), ...
    'Full and blocked periodic k-space apply paths should agree.');

io.assert_atomic_units(sys);
end

function sys = local_make_periodic_polsys()
sys = struct();

% Three active polarizable sites plus one nonpolarizable spectator. The
% mildly triclinic row-lattice prevents hidden row/column convention bugs.
sys.site_pos = [
     3.0   2.0   2.0
    10.0   6.0   5.0
    17.0  11.0   9.0
     6.0  15.0  12.0
];

sys.site_charge = zeros(4, 1);

% Small polarizabilities keep the fixed-point problem comfortably stable.
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

function local_assert_mu_close(mu, ref, tol, msg)
err = norm(mu - ref, 'fro');
scale = max(1, norm(ref, 'fro'));

assert(err <= tol * scale, ...
    '%s Relative Frobenius error = %.6e, tolerance = %.6e.', ...
    msg, err / scale, tol);
end