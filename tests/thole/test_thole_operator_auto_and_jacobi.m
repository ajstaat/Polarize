function test_thole_operator_auto_and_jacobi()
%TEST_THOLE_OPERATOR_AUTO_AND_JACOBI Verify operator auto resolution and Jacobi.
%
% Checks:
%   direct + auto + finite rcut -> dense_matrix / nonperiodic_paircache_dense
%   jacobi + auto + finite rcut -> matrix_free / nonperiodic_paircache_apply
%   jacobi + dense override     -> dense_matrix / nonperiodic_paircache_dense
%   matrix-free apply matches dense apply
%   Jacobi matrix-free matches direct dense
%   jacobi + auto + Rcut=Inf errors clearly
%   pair_cache public backend errors with rename guidance

polsys = local_make_four_site_polsys();

Eext = [
    0.010  0.000 0.000
    0.005  0.002 0.000
   -0.003  0.001 0.000
    0.002 -0.001 0.000
];

scfParams = struct();
scfParams.use_thole = true;
scfParams.softening = 0.0;
scfParams.rcut = 3.5;
scfParams.tol = 1e-11;
scfParams.maxIter = 500;
scfParams.mixing = 0.6;
scfParams.omega = 1.0;
scfParams.verbose = false;

problem = thole.prepare_scf_problem(polsys, Eext, scfParams);

opDirectAuto = thole.make_polarization_operator(polsys, problem, ...
    'Mode', 'nonperiodic', ...
    'Solver', 'direct', ...
    'Backend', 'auto', ...
    'UseThole', true, ...
    'Softening', 0.0, ...
    'Rcut', scfParams.rcut, ...
    'Verbose', false);

assert(strcmp(opDirectAuto.kind, 'dense_matrix'), ...
    'direct + auto should produce a dense_matrix operator.');

assert(strcmp(opDirectAuto.backend, 'nonperiodic_paircache_dense'), ...
    'direct + auto + finite rcut should use nonperiodic_paircache_dense.');

opJacobiAuto = thole.make_polarization_operator(polsys, problem, ...
    'Mode', 'nonperiodic', ...
    'Solver', 'jacobi', ...
    'Backend', 'auto', ...
    'UseThole', true, ...
    'Softening', 0.0, ...
    'Rcut', scfParams.rcut, ...
    'Verbose', false);

assert(strcmp(opJacobiAuto.kind, 'matrix_free'), ...
    'jacobi + auto should produce a matrix_free operator.');

assert(strcmp(opJacobiAuto.backend, 'nonperiodic_paircache_apply'), ...
    'jacobi + auto + finite rcut should use nonperiodic_paircache_apply.');

opJacobiDense = thole.make_polarization_operator(polsys, problem, ...
    'Mode', 'nonperiodic', ...
    'Solver', 'jacobi', ...
    'Backend', 'dense', ...
    'UseThole', true, ...
    'Softening', 0.0, ...
    'Rcut', scfParams.rcut, ...
    'Verbose', false);

assert(strcmp(opJacobiDense.kind, 'dense_matrix'), ...
    'jacobi + dense should produce a dense_matrix operator.');

assert(strcmp(opJacobiDense.backend, 'nonperiodic_paircache_dense'), ...
    'jacobi + dense + finite rcut should use nonperiodic_paircache_dense.');

muTest = (1:size(opDirectAuto.Tpol, 2)).' / 100;

assert(norm(opJacobiAuto.apply(muTest) - opDirectAuto.Tpol * muTest) < 1e-12, ...
    'matrix-free pair-cache apply should match dense Tpol action.');

assert(norm(opJacobiDense.apply(muTest) - opDirectAuto.Tpol * muTest) < 1e-12, ...
    'dense override apply should match direct dense Tpol action.');

[muDirect, infoDirect] = thole.solve_scf_direct(problem, opDirectAuto);

iterOpts = struct();
iterOpts.tol = 1e-11;
iterOpts.max_iter = 500;
iterOpts.mixing = 0.6;
iterOpts.verbose = false;

[muJacobi, infoJacobi] = thole.solve_scf_jacobi(problem, opJacobiAuto, iterOpts);

assert(infoDirect.relres < 1e-12, ...
    'Direct dense reference should have small residual.');

assert(infoJacobi.converged, ...
    'Jacobi matrix-free solver should converge.');

assert(infoJacobi.relres < 1e-11, ...
    'Jacobi matrix-free residual should be small.');

assert(norm(muJacobi - muDirect, 'fro') < 1e-9, ...
    'Jacobi matrix-free solution should match direct dense solution.');

energyDirect = calc.compute_total_energy_active_space(polsys, problem, muDirect, Eext, opDirectAuto);
energyJacobi = calc.compute_total_energy_active_space(polsys, problem, muJacobi, Eext, opDirectAuto);

assert(abs(energyJacobi.total - energyDirect.total) < 1e-10, ...
    'Jacobi energy should match direct energy.');

% Jacobi + auto + Rcut=Inf should error rather than silently building dense.
didError = false;

try
    thole.make_polarization_operator(polsys, problem, ...
        'Mode', 'nonperiodic', ...
        'Solver', 'jacobi', ...
        'Backend', 'auto', ...
        'Rcut', Inf);
catch ME
    didError = true;
    assert(strcmp(ME.identifier, ...
        'thole:make_polarization_operator:AutoMatrixFreeRequiresFiniteCutoff'), ...
        'jacobi + auto + Rcut=Inf should throw the finite-cutoff error.');
end

assert(didError, ...
    'jacobi + auto + Rcut=Inf should error.');

% Old pair_cache public backend should error with rename guidance.
didError = false;

try
    thole.make_polarization_operator(polsys, problem, ...
        'Mode', 'nonperiodic', ...
        'Solver', 'jacobi', ...
        'Backend', 'pair_cache', ...
        'Rcut', scfParams.rcut);
catch ME
    didError = true;
    assert(strcmp(ME.identifier, ...
        'thole:make_polarization_operator:PairCacheBackendRenamed'), ...
        'Backend="pair_cache" should throw the renamed-backend identifier.');
end

assert(didError, ...
    'Backend="pair_cache" should error with rename guidance.');

io.assert_atomic_units(polsys);

end

function polsys = local_make_four_site_polsys()

polsys = struct();

polsys.site_pos = [
    0.0 0.0 0.0
    3.0 0.0 0.0
    0.0 3.0 0.0
    6.0 0.0 0.0
];

polsys.site_charge = [0.0; 0.0; 0.0; 0.0];

% Small enough alphas for stable simple Jacobi convergence.
polsys.site_alpha = [0.20; 0.25; 0.22; 0.18];

polsys.site_is_polarizable = [true; true; true; true];

polsys.site_type = {'X'; 'X'; 'X'; 'X'};
polsys.site_class = {'p1'; 'p2'; 'p3'; 'p4'};
polsys.site_label = {'p1'; 'p2'; 'p3'; 'p4'};
polsys.site_mol_id = [1; 2; 3; 4];
polsys.site_is_active = [true; true; true; true];

polsys.n_sites = 4;

polsys.units = struct();
polsys.units.length = 'bohr';
polsys.units.alpha = 'atomic_unit';
polsys.units.charge = 'elementary_charge';

polsys.thole_a = 0.39;

polsys.is_periodic = false;
polsys.periodic_mode = 'nonperiodic';

end