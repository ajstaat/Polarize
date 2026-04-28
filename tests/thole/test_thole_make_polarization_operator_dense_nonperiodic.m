function test_thole_make_polarization_operator_dense_nonperiodic()
%TEST_THOLE_MAKE_POLARIZATION_OPERATOR_DENSE_NONPERIODIC Verify Operator-1.
%
% Checks:
%   - make_polarization_operator builds dense nonperiodic operator
%   - op.Tpol matches direct assembler
%   - op.apply matches Tpol * muVec
%   - solve_scf_direct accepts op or raw Tpol
%   - energy accepts op or raw Tpol
%   - unsupported backend errors clearly

polsys = local_make_two_site_polsys();

Eext = [
    0.10 0.0 0.0
    0.05 0.0 0.0
];

scfParams = struct();
scfParams.use_thole = true;
scfParams.softening = 0.0;
scfParams.rcut = Inf;
scfParams.tol = 1e-12;
scfParams.maxIter = 50;
scfParams.mixing = 1.0;
scfParams.omega = 1.0;
scfParams.verbose = false;

problem = thole.prepare_scf_problem(polsys, Eext, scfParams);

[Tref, infoRef] = thole.assemble_nonperiodic_interaction_matrix( ...
    polsys, problem, scfParams);

op = thole.make_polarization_operator(polsys, problem, ...
    'Mode', 'nonperiodic', ...
    'Backend', 'dense', ...
    'UseThole', scfParams.use_thole, ...
    'Softening', scfParams.softening, ...
    'Rcut', scfParams.rcut, ...
    'Verbose', false);

assert(strcmp(op.mode, 'nonperiodic'), ...
    'Operator mode should be nonperiodic.');

assert(strcmp(op.backend, 'dense'), ...
    'Operator backend should be dense.');

assert(strcmp(op.kind, 'dense_matrix'), ...
    'Dense backend should set op.kind = dense_matrix.');

assert(op.nPolSites == problem.nPolSites, ...
    'Operator nPolSites should match problem.');

assert(isequal(op.size, size(Tref)), ...
    'Operator size should match assembled Tpol size.');

assert(isfield(op, 'Tpol') && isnumeric(op.Tpol), ...
    'Dense operator should contain numeric op.Tpol.');

assert(norm(op.Tpol - Tref, 'fro') < 1e-12, ...
    'Operator Tpol should match direct nonperiodic assembler.');

assert(op.info.nPairBlocksKept == infoRef.nPairBlocksKept, ...
    'Operator info should preserve assembler diagnostics.');

muTest = (1:size(op.Tpol, 2)).' / 100;

assert(norm(op.apply(muTest) - op.Tpol * muTest) < 1e-12, ...
    'op.apply(muVec) should match op.Tpol * muVec.');

[muMatrix, infoMatrix] = thole.solve_scf_direct(problem, Tref);
[muOp, infoOp] = thole.solve_scf_direct(problem, op);

assert(norm(muMatrix - muOp, 'fro') < 1e-12, ...
    'Direct solver with op should match direct solver with raw Tpol.');

assert(abs(infoMatrix.relres - infoOp.relres) < 1e-14, ...
    'Direct solver residual should be unchanged by operator wrapper.');

energyMatrix = calc.compute_total_energy_active_space(polsys, problem, muMatrix, Eext, Tref);
energyOp = calc.compute_total_energy_active_space(polsys, problem, muOp, Eext, op);

assert(abs(energyMatrix.total - energyOp.total) < 1e-12, ...
    'Energy with op should match energy with raw Tpol.');

assert(abs(energyOp.stationary_consistency) < 1e-11, ...
    'Energy stationary consistency should remain small.');

% Unsupported backend should fail clearly.
didError = false;

try
    thole.make_polarization_operator(polsys, problem, ...
        'Mode', 'nonperiodic', ...
        'Backend', 'pair_cache');
catch ME
    didError = true;
    assert(strcmp(ME.identifier, ...
        'thole:make_polarization_operator:UnsupportedOperatorBackend'), ...
        'Unsupported backend should throw the expected identifier.');
end

assert(didError, ...
    'Unsupported pair_cache backend should error in Operator-1.');

io.assert_atomic_units(polsys);

end

function polsys = local_make_two_site_polsys()

polsys = struct();

polsys.site_pos = [
    0.0 0.0 0.0
    3.0 0.0 0.0
];

polsys.site_charge = [0.0; 0.0];

polsys.site_alpha = [1.0; 1.5];

polsys.site_is_polarizable = [true; true];

polsys.site_type = {'X'; 'X'};
polsys.site_class = {'p1'; 'p2'};
polsys.site_label = {'p1'; 'p2'};
polsys.site_mol_id = [1; 2];
polsys.site_is_active = [true; true];

polsys.n_sites = 2;

polsys.units = struct();
polsys.units.length = 'bohr';
polsys.units.alpha = 'atomic_unit';
polsys.units.charge = 'elementary_charge';

polsys.thole_a = 0.39;

polsys.is_periodic = false;
polsys.periodic_mode = 'nonperiodic';

end