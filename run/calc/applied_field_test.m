%% periodic_sor_toy_regression_test
% Regression test for periodic matrix-free SOR-like SCF on a tiny toy system.
%
% This compares:
%   1) direct solve using assembled periodic Tpol
%   2) matrix-free periodic SOR-like solve using
%        thole.solve_scf_iterative_periodic_sor(...)
%
% Goal:
%   Verify that the periodic matrix-free SOR solver reproduces the
%   assembled periodic reference solution on a tiny periodic test case.

clear; clc; close all;

fprintf('=== periodic toy SOR solver regression test ===\n');

%% ------------------------------------------------------------------------
% User controls

cfg = struct();

cfg.softening = 0.0;
cfg.use_thole = true;

cfg.ewald_tol = 1e-4;
cfg.boundary = 'tinfoil';
cfg.rcut_fraction = 0.7;

cfg.verbose = true;

cfg.omega = 0.9;
cfg.maxIter = 200;
cfg.tol = 1e-9;
cfg.printEvery = 10;
cfg.residualEvery = 5;
cfg.stopMetric = 'max_dmu';

% Tiny cubic toy cell
cfg.L_ang = 8.0;

% 4-site toy basis in Angstrom (converted below to bohr)
cfg.site_pos_ang = [ ...
    1.0, 1.0, 1.0; ...
    5.0, 1.0, 1.0; ...
    1.0, 5.0, 1.0; ...
    1.0, 1.0, 5.0];

cfg.site_alpha = [1.20; 1.10; 1.00; 0.90];
cfg.site_is_polarizable = true(4,1);
cfg.thole_a = 0.39;

tol = struct();
tol.direct_relres = 1e-12;
tol.sor_relres    = 1e-6;
tol.mu            = 1e-8;
tol.eHa           = 1e-10;

%% ------------------------------------------------------------------------
% Build tiny toy periodic polarization system directly in canonical atomic units

ANG2BOHR = 1.8897261254578281;

polsys = struct();
polsys.n_sites = size(cfg.site_pos_ang, 1);
polsys.site_pos = cfg.site_pos_ang * ANG2BOHR;        % bohr
polsys.site_alpha = cfg.site_alpha;                   % atomic units
polsys.site_is_polarizable = cfg.site_is_polarizable;
polsys.thole_a = cfg.thole_a;
polsys.site_charge = zeros(polsys.n_sites, 1);        % elementary charge
polsys.lattice = (cfg.L_ang * ANG2BOHR) * eye(3);     % bohr

polsys.units = struct();
polsys.units.length = 'bohr';
polsys.units.alpha = 'atomic_unit';
polsys.units.charge = 'elementary_charge';

io.assert_atomic_units(polsys);

fprintf('Toy periodic polarization system:\n');
fprintf('  n_sites             = %d\n', polsys.n_sites);
fprintf('  n_polarizable_sites = %d\n', nnz(polsys.site_is_polarizable));
fprintf('  box length          = %.6f bohr\n', cfg.L_ang * ANG2BOHR);

%% ------------------------------------------------------------------------
% External field
%
% Use a small nonzero external field so the SCF solve is nontrivial.

Eext = zeros(polsys.n_sites, 3);
Eext(:, 1) = 1.0e-3;
Eext(:, 2) = -5.0e-4;
Eext(:, 3) = 2.0e-4;

%% ------------------------------------------------------------------------
% SCF and Ewald parameters

scfParams = struct();
scfParams.softening = cfg.softening;
scfParams.use_thole = cfg.use_thole;
scfParams.verbose = cfg.verbose;
scfParams.printEvery = cfg.printEvery;
scfParams.residualEvery = cfg.residualEvery;
scfParams.omega = cfg.omega;
scfParams.maxIter = cfg.maxIter;
scfParams.tol = cfg.tol;
scfParams.stopMetric = cfg.stopMetric;

problem = thole.prepare_scf_problem(polsys, Eext, scfParams);

H = local_get_direct_lattice(polsys);
ewaldParams = ewald.choose_ewald_params(H, cfg.ewald_tol, cfg.boundary, cfg.rcut_fraction);

fprintf('\nChosen Ewald parameters:\n');
fprintf('  alpha    = %.6g\n', ewaldParams.alpha);
fprintf('  rcut     = %.6g\n', ewaldParams.rcut);
fprintf('  kcut     = %.6g\n', ewaldParams.kcut);
fprintf('  boundary = %s\n', ewaldParams.boundary);

%% ------------------------------------------------------------------------
% Assemble periodic operator

fprintf('\n--- assembled periodic operator ---\n');
tic;
[Tpol, parts, opinfo] = ewald.assemble_periodic_interaction_matrix( ...
    polsys, problem, ewaldParams, scfParams);
time_asm = toc;

fprintf('assembly time        = %.6f s\n', time_asm);
fprintf('operator size        = %d x %d\n', size(Tpol,1), size(Tpol,2));
fprintf('nRealInteractions    = %d\n', opinfo.nRealInteractions);
fprintf('num_kvec             = %d\n', opinfo.num_kvec);

%% ------------------------------------------------------------------------
% Direct reference solve

fprintf('\n--- direct periodic solve ---\n');
tic;
[mu_direct, direct] = thole.solve_scf_direct(problem, Tpol);
time_direct = toc;

E_direct = calc.compute_total_energy_active_space(polsys, problem, mu_direct, Eext, Tpol);
relres_direct = thole.compute_active_space_relres(problem, Tpol, mu_direct);

fprintf('direct: relres = %.3e | E = %+0.12f Ha | time = %.6f s\n', ...
    relres_direct, E_direct.total, time_direct);

%% ------------------------------------------------------------------------
% Matrix-free periodic SOR solve

fprintf('\n--- matrix-free periodic SOR solve ---\n');
tic;
[mu_sor, scf_sor] = thole.solve_scf_iterative_periodic_sor( ...
    polsys, Eext, ewaldParams, scfParams);
time_sor = toc;

problemSOR = thole.prepare_scf_problem(polsys, Eext, scfParams);
E_sor = calc.compute_total_energy_active_space(polsys, problemSOR, mu_sor, Eext, Tpol);
relres_sor = thole.compute_active_space_relres(problemSOR, Tpol, mu_sor);

fprintf('iter-SOR: relres = %.3e | E = %+0.12f Ha | iters = %d | time = %.6f s | converged = %d\n', ...
    relres_sor, E_sor.total, scf_sor.nIter, time_sor, scf_sor.converged);

%% ------------------------------------------------------------------------
% Comparisons

fprintf('\n--- comparisons ---\n');

mu_diff = max(sqrt(sum((mu_direct - mu_sor).^2, 2)));
e_diff  = abs(E_direct.total - E_sor.total);

fprintf('max |mu_direct - mu_sor| = %.16e\n', mu_diff);
fprintf('|E_direct - E_sor|       = %.16e Ha\n', e_diff);

%% ------------------------------------------------------------------------
% Assertions

if ~isfinite(relres_direct) || relres_direct > tol.direct_relres
    error('Direct periodic solver relative residual too large: %.3e', relres_direct);
end

if ~scf_sor.converged
    error('Periodic matrix-free SOR solver did not converge.');
end

if ~isfinite(relres_sor) || relres_sor > tol.sor_relres
    error('Periodic matrix-free SOR relative residual too large: %.3e', relres_sor);
end

if mu_diff > tol.mu
    error('Direct and periodic matrix-free SOR dipoles disagree: %.3e', mu_diff);
end

if e_diff > tol.eHa
    error('Direct and periodic matrix-free SOR energies disagree: %.3e Ha', e_diff);
end

fprintf('\nPASS: periodic matrix-free SOR solver matches direct periodic reference.\n');

%% ========================================================================
% Local helpers
%% ========================================================================

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('local_get_direct_lattice:MissingLattice', ...
            'Need sys.super_lattice or sys.lattice.');
    end
end