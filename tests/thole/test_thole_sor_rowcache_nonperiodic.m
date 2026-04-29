function test_thole_sor_rowcache_nonperiodic()
%TEST_THOLE_SOR_ROWCACHE_NONPERIODIC Verify SOR with row-cache operator.

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

opDense = thole.make_polarization_operator(polsys, problem, ...
    'Mode', 'nonperiodic', ...
    'Solver', 'direct', ...
    'Backend', 'auto', ...
    'UseThole', true, ...
    'Softening', 0.0, ...
    'Rcut', scfParams.rcut, ...
    'Verbose', false);

opSor = thole.make_polarization_operator(polsys, problem, ...
    'Mode', 'nonperiodic', ...
    'Solver', 'sor', ...
    'Backend', 'auto', ...
    'UseThole', true, ...
    'Softening', 0.0, ...
    'Rcut', scfParams.rcut, ...
    'Verbose', false);

assert(strcmp(opSor.kind, 'matrix_free'), ...
    'SOR auto operator should be matrix_free.');

assert(strcmp(opSor.backend, 'nonperiodic_rowcache_apply'), ...
    'SOR auto operator should use nonperiodic_rowcache_apply backend.');

assert(opSor.capabilities.apply, ...
    'SOR operator should support full apply.');

assert(opSor.capabilities.row_update, ...
    'SOR operator should support row_update.');

assert(~opSor.capabilities.dense_matrix, ...
    'SOR row-cache operator should not report dense_matrix capability.');

assert(isfield(opSor, 'apply_row') && isa(opSor.apply_row, 'function_handle'), ...
    'SOR row-cache operator should provide op.apply_row.');

muTest = (1:size(opDense.Tpol, 2)).' / 100;

assert(norm(opSor.apply(muTest) - opDense.Tpol * muTest) < 1e-12, ...
    'SOR row-cache full apply should match dense Tpol action.');

for rowLocal = 1:problem.nPolSites
    block = (3*(rowLocal-1)+1):(3*rowLocal);
    rowDense = opDense.Tpol(block, :) * muTest;
    rowCache = opSor.apply_row(rowLocal, muTest);

    assert(norm(rowCache - rowDense) < 1e-12, ...
        'SOR row-cache apply_row should match dense row action.');
end

[muDirect, infoDirect] = thole.solve_scf_direct(problem, opDense);

sorOpts = struct();
sorOpts.tol = 1e-11;
sorOpts.max_iter = 500;
sorOpts.omega = 1.0;
sorOpts.verbose = false;

[muSor, infoSor] = thole.solve_scf_sor(problem, opSor, sorOpts);

assert(infoDirect.relres < 1e-12, ...
    'Direct dense reference should have small residual.');

assert(infoSor.converged, ...
    'SOR row-cache solver should converge.');

assert(infoSor.relres < 1e-11, ...
    'SOR row-cache residual should be small.');

assert(norm(muSor - muDirect, 'fro') < 1e-9, ...
    'SOR row-cache solution should match dense direct solution.');

energyDirect = calc.compute_total_energy_active_space(polsys, problem, muDirect, Eext, opDense);
energySor = calc.compute_total_energy_active_space(polsys, problem, muSor, Eext, opDense);

assert(abs(energySor.total - energyDirect.total) < 1e-10, ...
    'SOR energy should match direct energy.');

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