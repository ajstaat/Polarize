function test_thole_gmres_matrix_free_nonperiodic()
%TEST_THOLE_GMRES_MATRIX_FREE_NONPERIODIC Verify GMRES with matrix-free op.
%
% Checks:
%   - gmres + auto + finite rcut builds matrix-free pair-cache apply op
%   - matrix-free apply matches dense apply
%   - GMRES(matrix_free) matches Direct(dense)
%   - energy from GMRES solution matches direct energy

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
scfParams.tol = 1e-12;
scfParams.maxIter = 100;
scfParams.mixing = 0.6;
scfParams.omega = 1.0;
scfParams.verbose = false;

problem = thole.prepare_scf_problem(polsys, Eext, scfParams);

opDense = thole.make_polarization_operator(polsys, problem, ...
    'Mode', 'nonperiodic', ...
    'Solver', 'direct', ...
    'Backend', 'auto', ...
    'UseThole', scfParams.use_thole, ...
    'Softening', scfParams.softening, ...
    'Rcut', scfParams.rcut, ...
    'Verbose', false);

assert(strcmp(opDense.kind, 'dense_matrix'), ...
    'direct + auto should build a dense matrix operator.');

assert(strcmp(opDense.backend, 'nonperiodic_paircache_dense'), ...
    'direct + auto finite-cutoff backend should be nonperiodic_paircache_dense.');

opGmres = thole.make_polarization_operator(polsys, problem, ...
    'Mode', 'nonperiodic', ...
    'Solver', 'gmres', ...
    'Backend', 'auto', ...
    'UseThole', scfParams.use_thole, ...
    'Softening', scfParams.softening, ...
    'Rcut', scfParams.rcut, ...
    'Verbose', false);

assert(strcmp(opGmres.kind, 'matrix_free'), ...
    'gmres + auto should build a matrix-free operator.');

assert(strcmp(opGmres.backend, 'nonperiodic_paircache_apply'), ...
    'gmres + auto finite-cutoff backend should be nonperiodic_paircache_apply.');

muTest = (1:size(opDense.Tpol, 2)).' / 100;

assert(norm(opGmres.apply(muTest) - opDense.Tpol * muTest) < 1e-12, ...
    'GMRES matrix-free operator action should match dense Tpol action.');

[muDirect, infoDirect] = thole.solve_scf_direct(problem, opDense);

gmresOpts = struct();
gmresOpts.tol = 1e-12;
gmresOpts.max_iter = 50;
gmresOpts.restart = [];
gmresOpts.verbose = false;

[muGmres, infoGmres] = thole.solve_scf_gmres(problem, opGmres, gmresOpts);

assert(infoDirect.relres < 1e-12, ...
    'Direct dense reference should have small residual.');

assert(infoGmres.converged, ...
    'GMRES matrix-free solver should converge.');

assert(infoGmres.relres < 1e-11, ...
    'GMRES matrix-free residual should be small.');

assert(norm(muGmres - muDirect, 'fro') < 1e-9, ...
    'GMRES matrix-free solution should match dense direct solution.');

energyDirect = calc.compute_total_energy_active_space(polsys, problem, muDirect, Eext, opDense);
energyGmres = calc.compute_total_energy_active_space(polsys, problem, muGmres, Eext, opDense);

assert(abs(energyGmres.total - energyDirect.total) < 1e-10, ...
    'GMRES energy should match direct energy.');

% Explicit dense backend should also work with GMRES.
opGmresDense = thole.make_polarization_operator(polsys, problem, ...
    'Mode', 'nonperiodic', ...
    'Solver', 'gmres', ...
    'Backend', 'dense', ...
    'UseThole', scfParams.use_thole, ...
    'Softening', scfParams.softening, ...
    'Rcut', scfParams.rcut, ...
    'Verbose', false);

assert(strcmp(opGmresDense.kind, 'dense_matrix'), ...
    'gmres + dense should build a dense matrix operator.');

[muGmresDense, infoGmresDense] = thole.solve_scf_gmres(problem, opGmresDense, gmresOpts);

assert(infoGmresDense.converged, ...
    'GMRES dense operator solver should converge.');

assert(norm(muGmresDense - muDirect, 'fro') < 1e-9, ...
    'GMRES dense operator solution should match direct solution.');

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

% Small/stable alphas.
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