function op = build_periodic_dense_operator(sys, problem, opt)
%BUILD_PERIODIC_DENSE_OPERATOR Build dense periodic Ewald polarization op.
%
% Private backend builder used by thole.make_polarization_operator.
%
% Returns:
%   op.kind    = 'dense_matrix'
%   op.backend = 'periodic_paircache_dense'
%
% This dense backend is intended for tiny validation systems and direct SCF.
% It constructs the dense active-space matrix by applying the periodic
% matrix-free paircache operator to Cartesian basis vectors.
%
% For production periodic calculations, use:
%   - build_periodic_paircache_operator for Jacobi/GMRES
%   - build_periodic_rowcache_operator for SOR

io.assert_atomic_units(sys);

nPol = problem.nPolSites;
nDim = 3 * nPol;

if nDim == 0
    Tpol = zeros(0, 0);
    pairOp = build_periodic_paircache_operator(sys, problem, opt);

    op = struct();
    op.mode = 'periodic_ewald';
    op.kind = 'dense_matrix';
    op.backend = 'periodic_paircache_dense';
    op.nPolSites = nPol;
    op.size = [0, 0];
    op.Tpol = Tpol;
    op.apply = @(muVec) local_apply_dense(muVec, Tpol, nPol);
    op.info = pairOp.info;
    op.info.assembly_backend = 'periodic_paircache_dense';
    op.info.dense_build_time = 0.0;
    op.info.dense_source_backend = pairOp.backend;

    op.capabilities = struct();
    op.capabilities.apply = true;
    op.capabilities.dense_matrix = true;
    op.capabilities.row_update = false;

    op.params = pairOp.params;
    op.request = struct();

    return;
end

pairOp = build_periodic_paircache_operator(sys, problem, opt);

tDense = tic;

Tpol = zeros(nDim, nDim);

for j = 1:nDim
    ej = zeros(nDim, 1);
    ej(j) = 1.0;
    Tpol(:, j) = pairOp.apply(ej);
end

denseBuildTime = toc(tDense);

op = struct();

op.mode = 'periodic_ewald';
op.kind = 'dense_matrix';
op.backend = 'periodic_paircache_dense';

op.nPolSites = nPol;
op.size = size(Tpol);

op.Tpol = Tpol;
op.apply = @(muVec) local_apply_dense(muVec, Tpol, nPol);

op.info = pairOp.info;
op.info.assembly_backend = 'periodic_paircache_dense';
op.info.dense_build_time = denseBuildTime;
op.info.dense_source_backend = pairOp.backend;
op.info.dense_nDim = nDim;

op.cache = pairOp.cache;
op.real_cache = pairOp.real_cache;
op.k_cache = pairOp.k_cache;
op.periodic_cache = pairOp.periodic_cache;

op.capabilities = struct();
op.capabilities.apply = true;
op.capabilities.dense_matrix = true;
op.capabilities.row_update = false;

op.params = pairOp.params;
end

function y = local_apply_dense(muVec, Tpol, nPol)
muVec = muVec(:);

if numel(muVec) ~= 3*nPol
    error('thole:build_periodic_dense_operator:BadApplyVectorSize', ...
        'muVec must have length 3*nPolSites.');
end

y = Tpol * muVec;
end