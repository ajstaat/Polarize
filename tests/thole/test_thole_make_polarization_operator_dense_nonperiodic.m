function test_thole_make_polarization_operator_dense_nonperiodic()
%TEST_THOLE_MAKE_POLARIZATION_OPERATOR_DENSE_NONPERIODIC Verify dense Operator-1.
%
% Checks:
%   - make_polarization_operator builds dense nonperiodic operator
%   - finite Rcut dense operator uses nonperiodic_paircache_dense backend
%   - op.Tpol matches a brute-force finite-cutoff reference
%   - op.apply matches op.Tpol * muVec
%   - solve_scf_direct accepts op or raw Tpol
%   - energy accepts op or raw Tpol
%   - rcut = Inf still uses all-pairs dense backend
%   - unsupported/renamed backend errors clearly

polsys = local_make_four_site_polsys();

Eext = [
    0.10  0.00 0.00
    0.05  0.02 0.00
   -0.03  0.01 0.00
    0.02 -0.01 0.00
];

scfParams = struct();
scfParams.use_thole = true;
scfParams.softening = 0.0;
scfParams.rcut = 3.5;
scfParams.use_mex = true;
scfParams.profile = false;
scfParams.tol = 1e-12;
scfParams.maxIter = 50;
scfParams.mixing = 1.0;
scfParams.omega = 1.0;
scfParams.verbose = false;

problem = thole.prepare_scf_problem(polsys, Eext, scfParams);

% -------------------------------------------------------------------------
% Finite-cutoff dense operator should use pair-cache dense backend.
% -------------------------------------------------------------------------

op = thole.make_polarization_operator(polsys, problem, ...
    'Mode', 'nonperiodic', ...
    'Solver', 'direct', ...
    'Backend', 'dense', ...
    'UseThole', scfParams.use_thole, ...
    'Softening', scfParams.softening, ...
    'Rcut', scfParams.rcut, ...
    'UseMex', scfParams.use_mex, ...
    'Profile', scfParams.profile, ...
    'Verbose', false);

assert(strcmp(op.mode, 'nonperiodic'), ...
    'Operator mode should be nonperiodic.');

assert(strcmp(op.kind, 'dense_matrix'), ...
    'Dense operator should set op.kind = dense_matrix.');

assert(strcmp(op.backend, 'nonperiodic_paircache_dense'), ...
    'Finite-cutoff dense operator should use nonperiodic_paircache_dense backend.');

assert(op.nPolSites == problem.nPolSites, ...
    'Operator nPolSites should match problem.');

assert(isequal(op.size, [3*problem.nPolSites 3*problem.nPolSites]), ...
    'Operator size should be 3*nPolSites square.');

assert(isfield(op, 'Tpol') && isnumeric(op.Tpol), ...
    'Dense operator should contain numeric op.Tpol.');

assert(isfield(op, 'apply') && isa(op.apply, 'function_handle'), ...
    'Dense operator should provide op.apply.');

assert(isfield(op, 'info'), ...
    'Operator should contain assembly diagnostics.');

assert(strcmp(op.info.assembly_backend, 'pair_cache'), ...
    'Finite cutoff dense assembly should use pair_cache assembly backend.');

assert(op.info.use_cutoff, ...
    'Finite cutoff operator should report use_cutoff = true.');

assert(op.info.rcut == scfParams.rcut, ...
    'Operator info should preserve rcut.');

assert(op.info.nPairBlocksKept == 3, ...
    'For this geometry and rcut, exactly three polarizable pairs should be kept.');

assert(op.info.nPairBlocksSkippedCutoff == 3, ...
    'For four polarizable sites, three of six pairs should be skipped by cutoff.');

assert(isfield(op.info, 'cache_info'), ...
    'Pair-cache dense assembly should include cache_info diagnostics.');

assert(op.info.cache_info.n_pairs_returned == op.info.nPairBlocksKept, ...
    'Cache pair count should match kept pair-block count.');

assert(isfield(op, 'capabilities'), ...
    'Operator should contain capabilities.');

assert(op.capabilities.apply, ...
    'Dense operator should report apply capability.');

assert(op.capabilities.dense_matrix, ...
    'Dense operator should report dense_matrix capability.');

assert(~op.capabilities.row_update, ...
    'Dense operator should not report row_update capability.');

% -------------------------------------------------------------------------
% Compare against brute-force finite-cutoff reference.
% -------------------------------------------------------------------------

Tref = local_build_bruteforce_Tpol(polsys, problem, scfParams);

assert(norm(op.Tpol - Tref, 'fro') < 1e-12, ...
    'Pair-cache/vectorized finite-cutoff Tpol should match brute-force reference.');

muTest = (1:size(op.Tpol, 2)).' / 100;

assert(norm(op.apply(muTest) - op.Tpol * muTest) < 1e-12, ...
    'op.apply(muVec) should match op.Tpol * muVec.');

% -------------------------------------------------------------------------
% Direct solver and energy should accept either raw Tpol or operator.
% -------------------------------------------------------------------------

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

% -------------------------------------------------------------------------
% Infinite cutoff should use all-pairs dense backend.
% -------------------------------------------------------------------------

opInf = thole.make_polarization_operator(polsys, problem, ...
    'Mode', 'nonperiodic', ...
    'Solver', 'direct', ...
    'Backend', 'dense', ...
    'UseThole', scfParams.use_thole, ...
    'Softening', scfParams.softening, ...
    'Rcut', Inf, ...
    'UseMex', scfParams.use_mex, ...
    'Profile', scfParams.profile, ...
    'Verbose', false);

assert(strcmp(opInf.mode, 'nonperiodic'), ...
    'Infinite-cutoff operator mode should be nonperiodic.');

assert(strcmp(opInf.kind, 'dense_matrix'), ...
    'Infinite-cutoff operator should still be dense_matrix.');

assert(strcmp(opInf.backend, 'nonperiodic_allpairs_dense'), ...
    'Infinite-cutoff dense operator should use nonperiodic_allpairs_dense backend.');

assert(strcmp(opInf.info.assembly_backend, 'all_pairs'), ...
    'Infinite cutoff should use all-pairs assembly backend.');

assert(~opInf.info.use_cutoff, ...
    'Infinite cutoff operator should report use_cutoff = false.');

assert(opInf.info.nPairBlocksKept == 6, ...
    'Four polarizable sites should have six all-pairs blocks.');

assert(opInf.info.nPairBlocksSkippedCutoff == 0, ...
    'Infinite cutoff should skip no pair blocks.');

TrefInf = local_build_bruteforce_Tpol(polsys, problem, struct( ...
    'use_thole', true, ...
    'softening', 0.0, ...
    'rcut', Inf));

assert(norm(opInf.Tpol - TrefInf, 'fro') < 1e-12, ...
    'Infinite-cutoff operator should match all-pairs brute-force reference.');

% -------------------------------------------------------------------------
% Direct + auto should also resolve to dense.
% -------------------------------------------------------------------------

opAuto = thole.make_polarization_operator(polsys, problem, ...
    'Mode', 'nonperiodic', ...
    'Solver', 'direct', ...
    'Backend', 'auto', ...
    'UseThole', scfParams.use_thole, ...
    'Softening', scfParams.softening, ...
    'Rcut', scfParams.rcut, ...
    'Verbose', false);

assert(strcmp(opAuto.kind, 'dense_matrix'), ...
    'direct + auto should resolve to a dense_matrix operator.');

assert(strcmp(opAuto.backend, 'nonperiodic_paircache_dense'), ...
    'direct + auto + finite rcut should resolve to nonperiodic_paircache_dense.');

assert(norm(opAuto.Tpol - op.Tpol, 'fro') < 1e-12, ...
    'direct + auto dense operator should match explicit dense operator.');

% -------------------------------------------------------------------------
% Direct + matrix_free should fail clearly.
% -------------------------------------------------------------------------

didError = false;

try
    thole.make_polarization_operator(polsys, problem, ...
        'Mode', 'nonperiodic', ...
        'Solver', 'direct', ...
        'Backend', 'matrix_free', ...
        'Rcut', scfParams.rcut);
catch ME
    didError = true;
    assert(strcmp(ME.identifier, ...
        'thole:make_polarization_operator:DirectRequiresDenseOperator'), ...
        'direct + matrix_free should throw the expected identifier.');
end

assert(didError, ...
    'direct + matrix_free should error.');

% -------------------------------------------------------------------------
% Old public pair_cache backend name should fail with rename guidance.
% -------------------------------------------------------------------------

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

% =========================================================================
% Reference helpers
% =========================================================================

function Tpol = local_build_bruteforce_Tpol(sys, problem, scfParams)

sites = problem.activeSites(:);
nPol = problem.nPolSites;

Tpol = zeros(3*nPol, 3*nPol);

opts = struct();
opts.use_thole = scfParams.use_thole;
opts.softening = scfParams.softening;

if isfield(scfParams, 'rcut') && ~isempty(scfParams.rcut)
    rcut = scfParams.rcut;
else
    rcut = Inf;
end

for a = 1:(nPol - 1)
    i = sites(a);
    ri = sys.site_pos(i, :);

    Ia = local_block_indices(a);

    for b = (a + 1):nPol
        j = sites(b);
        rj = sys.site_pos(j, :);

        if isfinite(rcut)
            rij = ri - rj;
            if sum(rij.^2) > rcut^2
                continue;
            end
        end

        Ib = local_block_indices(b);

        Tij = thole.dipole_tensor_block( ...
            sys.site_pos(i, :), ...
            sys.site_pos(j, :), ...
            sys.site_alpha(i), ...
            sys.site_alpha(j), ...
            sys.thole_a, ...
            opts);

        Tpol(Ia, Ib) = Tij;
        Tpol(Ib, Ia) = Tij.';
    end
end

end

function idx = local_block_indices(k)

idx = (3*(k-1)+1):(3*k);

end

% =========================================================================
% Test system
% =========================================================================

function polsys = local_make_four_site_polsys()

polsys = struct();

% Four polarizable sites. With rcut = 3.5, the included pairs are:
%   1--2 distance 3
%   1--3 distance 3
%   2--4 distance 3
%
% The excluded pairs are:
%   1--4 distance 6
%   2--3 distance sqrt(18)
%   3--4 distance sqrt(45)
%
% This gives a useful finite-cutoff pattern that is not all-pairs.
polsys.site_pos = [
    0.0 0.0 0.0
    3.0 0.0 0.0
    0.0 3.0 0.0
    6.0 0.0 0.0
];

polsys.site_charge = [0.0; 0.0; 0.0; 0.0];

polsys.site_alpha = [1.0; 1.5; 1.2; 1.1];

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