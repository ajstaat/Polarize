%% tutorial_05_nonperiodic_solver_interface
% Walk through the nonperiodic polarization solver interface:
%
%   polsys
%   -> Eext
%   -> problem
%   -> operator
%   -> solver
%   -> energy
%
% This tutorial compares:
%
%   direct  : dense_matrix operator
%   jacobi  : matrix_free pair-cache operator
%   gmres   : matrix_free pair-cache operator
%   sor     : matrix_free row-cache operator
%
% The key design:
%
%   thole.make_polarization_operator(..., 'Solver', method, 'Backend', 'auto')
%
% chooses the appropriate operator representation for the solver.

clear; clc; close all;

fprintf('\n============================================================\n');
fprintf('Tutorial 05: Nonperiodic solver interface\n');
fprintf('============================================================\n');

HARTREE_TO_EV = 27.211386245988;

%% 1. Build a tiny nonperiodic polarization system

fprintf('\n[1] Building tiny nonperiodic polsys...\n');

polsys = local_make_demo_polsys();
io.assert_atomic_units(polsys);

fprintf('  n_sites             = %d\n', polsys.n_sites);
fprintf('  polarizable sites   = %d\n', nnz(polsys.site_is_polarizable));
fprintf('  charged sites       = %d\n', nnz(abs(polsys.site_charge) > 0));
fprintf('  length units        = %s\n', polsys.units.length);
fprintf('  alpha units         = %s\n', polsys.units.alpha);

%% 2. Compute nonperiodic external field from source charges

fprintf('\n[2] Computing nonperiodic external field...\n');

fieldParams = struct();
fieldParams.use_thole = false;
fieldParams.field = struct();
fieldParams.field.mode = 'nonperiodic';
fieldParams.field.exclude_self = true;
fieldParams.field.use_thole_damping = false;
fieldParams.field.target_mask = logical(polsys.site_is_polarizable(:));
fieldParams.field.source_mask = abs(polsys.site_charge(:)) > 0;

Eext = calc.compute_external_field(polsys, fieldParams);

fprintf('  ||Eext||_F              = %.12e\n', norm(Eext, 'fro'));
fprintf('  ||Eext polarizable||_F  = %.12e\n', ...
    norm(Eext(polsys.site_is_polarizable, :), 'fro'));

assert(norm(Eext(polsys.site_is_polarizable, :), 'fro') > 0, ...
    'Expected nonzero external field on polarizable sites.');

%% 3. Prepare active-space SCF problem

fprintf('\n[3] Preparing SCF problem...\n');

scfParams = struct();
scfParams.use_thole = true;
scfParams.softening = 0.0;
scfParams.rcut = 5.0;
scfParams.tol = 1e-11;
scfParams.maxIter = 500;
scfParams.mixing = 0.6;
scfParams.omega = 1.0;
scfParams.stop_metric = 'relres';
scfParams.verbose = false;

problem = thole.prepare_scf_problem(polsys, Eext, scfParams);

fprintf('  nPolSites          = %d\n', problem.nPolSites);
fprintf('  active vector size = %d\n', numel(problem.Eext_pol_vec));
fprintf('  rcut               = %.6g bohr\n', scfParams.rcut);

%% 4. Solve with each method

fprintf('\n[4] Solving with direct / Jacobi / GMRES / SOR...\n');

methods = {'direct', 'jacobi', 'gmres', 'sor'};
results = struct();

for k = 1:numel(methods)
    method = methods{k};

    fprintf('\n--- solver: %s ---\n', method);

    op = thole.make_polarization_operator(polsys, problem, ...
        'Mode', 'nonperiodic', ...
        'Solver', method, ...
        'Backend', 'auto', ...
        'UseThole', scfParams.use_thole, ...
        'Softening', scfParams.softening, ...
        'Rcut', scfParams.rcut, ...
        'UseMex', true, ...
        'Profile', false, ...
        'Verbose', false);

    fprintf('  op.kind       = %s\n', op.kind);
    fprintf('  op.backend    = %s\n', op.backend);

    tSolve = tic;

    switch method
        case 'direct'
            solveOpts = struct();
            solveOpts.compute_residual = true;

            [mu, info] = thole.solve_scf_direct(problem, op, solveOpts);

        case 'jacobi'
            solveOpts = struct();
            solveOpts.tol = scfParams.tol;
            solveOpts.max_iter = scfParams.maxIter;
            solveOpts.mixing = scfParams.mixing;
            solveOpts.stop_metric = scfParams.stop_metric;
            solveOpts.verbose = false;

            [mu, info] = thole.solve_scf_jacobi(problem, op, solveOpts);

        case 'gmres'
            solveOpts = struct();
            solveOpts.tol = scfParams.tol;
            solveOpts.max_iter = scfParams.maxIter;
            solveOpts.restart = [];
            solveOpts.verbose = false;

            [mu, info] = thole.solve_scf_gmres(problem, op, solveOpts);

        case 'sor'
            solveOpts = struct();
            solveOpts.tol = scfParams.tol;
            solveOpts.max_iter = scfParams.maxIter;
            solveOpts.omega = scfParams.omega;
            solveOpts.stop_metric = scfParams.stop_metric;
            solveOpts.verbose = false;

            [mu, info] = thole.solve_scf_sor(problem, op, solveOpts);

        otherwise
            error('Unknown method "%s".', method);
    end

    wallTime = toc(tSolve);

    % Energy currently requires dense Tpol, so use a dense reference op for
    % non-direct matrix-free solutions.
    if strcmp(op.kind, 'dense_matrix')
        energyOp = op;
    else
        energyOp = thole.make_polarization_operator(polsys, problem, ...
            'Mode', 'nonperiodic', ...
            'Solver', 'direct', ...
            'Backend', 'dense', ...
            'UseThole', scfParams.use_thole, ...
            'Softening', scfParams.softening, ...
            'Rcut', scfParams.rcut, ...
            'UseMex', true, ...
            'Profile', false, ...
            'Verbose', false);
    end

    energy = calc.compute_total_energy_active_space(polsys, problem, mu, Eext, energyOp);

    fprintf('  converged     = %d\n', logical(info.converged));
    fprintf('  relres        = %.12e\n', info.relres);
    fprintf('  solve wall    = %.6f s\n', wallTime);
    fprintf('  ||mu||_F      = %.12e\n', norm(mu, 'fro'));
    fprintf('  Epol          = %+ .12e Ha  (%+ .8f eV)\n', ...
        energy.total, energy.total * HARTREE_TO_EV);

    results.(method).op = op;
    results.(method).mu = mu;
    results.(method).info = info;
    results.(method).energy = energy;
    results.(method).wall_time = wallTime;
end

%% 5. Compare solver results against direct

fprintf('\n[5] Comparing against direct reference...\n');

muRef = results.direct.mu;
ERef = results.direct.energy.total;

fprintf('\n%-10s %14s %14s %14s %14s\n', ...
    'solver', 'relres', 'dmu_F', 'Epol/eV', 'dE/meV');
fprintf('%s\n', repmat('-', 1, 72));

for k = 1:numel(methods)
    method = methods{k};

    mu = results.(method).mu;
    energy = results.(method).energy;

    dmu = norm(mu - muRef, 'fro');
    dEmeV = (energy.total - ERef) * HARTREE_TO_EV * 1000;

    fprintf('%-10s %14.6e %14.6e %14.8f %14.6e\n', ...
        method, ...
        results.(method).info.relres, ...
        dmu, ...
        energy.total * HARTREE_TO_EV, ...
        dEmeV);
end

fprintf('\nTutorial completed successfully.\n');

%% Local helpers

function polsys = local_make_demo_polsys()

polsys = struct();

% Two charged source sites and four polarizable target sites.
polsys.site_pos = [
    0.0  0.0 0.0   % + charge source
    5.0  0.0 0.0   % - charge source
    1.5  2.0 0.0   % polarizable
    3.0  2.5 0.0   % polarizable
    4.5  2.0 0.0   % polarizable
    2.5 -2.0 0.0   % polarizable
];

polsys.site_charge = [
    +1.0
    -1.0
     0.0
     0.0
     0.0
     0.0
];

polsys.site_alpha = [
    0.0
    0.0
    0.40
    0.35
    0.42
    0.38
];

polsys.site_is_polarizable = [
    false
    false
    true
    true
    true
    true
];

polsys.site_type = {'X'; 'X'; 'X'; 'X'; 'X'; 'X'};
polsys.site_class = {'qplus'; 'qminus'; 'pol'; 'pol'; 'pol'; 'pol'};
polsys.site_label = {'q+'; 'q-'; 'p1'; 'p2'; 'p3'; 'p4'};
polsys.site_mol_id = (1:6).';
polsys.site_is_active = abs(polsys.site_charge) > 0;

polsys.n_sites = 6;

polsys.units = struct();
polsys.units.length = 'bohr';
polsys.units.alpha = 'atomic_unit';
polsys.units.charge = 'elementary_charge';

polsys.thole_a = 0.39;

polsys.is_periodic = false;
polsys.periodic_mode = 'nonperiodic';

end