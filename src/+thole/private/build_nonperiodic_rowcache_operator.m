function op = build_nonperiodic_rowcache_operator(sys, problem, opt)
%BUILD_NONPERIODIC_ROWCACHE_OPERATOR Build matrix-free row-cache operator.
%
% Private backend builder used by thole.make_polarization_operator.
%
% Returns:
%   op.kind    = 'matrix_free'
%   op.backend = 'nonperiodic_rowcache_apply'
%
% This backend stores directed row-cache tensor coefficients and exposes:
%   op.apply(muVec)
%   op.apply_row(rowLocal, muVec)
%
% It does not build or store dense op.Tpol.

if ~isfinite(opt.Rcut)
    error('thole:build_nonperiodic_rowcache_operator:RequiresFiniteCutoff', ...
        ['Nonperiodic row-cache matrix-free operator requires finite Rcut. ', ...
         'Use Backend="dense" for all-pairs dense calculations.']);
end

io.assert_atomic_units(sys);

nPol = problem.nPolSites;
nSites = size(sys.site_pos, 1);

cacheOpts = struct();
cacheOpts.rcut = opt.Rcut;
cacheOpts.use_mex = opt.UseMex;
cacheOpts.profile = opt.Profile || opt.Verbose;
cacheOpts.use_thole = opt.UseThole;

tCache = tic;
rowCache = geom.build_active_row_cache_nonperiodic(sys, problem, cacheOpts);
cacheTime = toc(tCache);

rowPtr = local_normalize_row_ptr(rowCache.row_ptr(:));
colIdxRaw = local_normalize_col_idx(rowCache.col_idx(:));

if numel(rowPtr) ~= nPol + 1
    error('thole:build_nonperiodic_rowcache_operator:BadRowPtr', ...
        'rowCache.row_ptr must have length nPolSites + 1.');
end

nEntries = rowPtr(end) - 1;

if numel(colIdxRaw) ~= nEntries
    error('thole:build_nonperiodic_rowcache_operator:BadColIdxLength', ...
        'numel(rowCache.col_idx) must equal rowPtr(end)-1.');
end

colIdx = local_map_columns_to_active(colIdxRaw, sys, problem);

if size(rowCache.dr, 1) ~= nEntries || size(rowCache.dr, 2) ~= 3
    error('thole:build_nonperiodic_rowcache_operator:BadDr', ...
        'rowCache.dr must be nEntries x 3.');
end

dr = double(rowCache.dr);

dx = dr(:, 1);
dy = dr(:, 2);
dz = dr(:, 3);

if isfield(rowCache, 'r2_bare') && ~isempty(rowCache.r2_bare)
    r2Bare = double(rowCache.r2_bare(:));
else
    r2Bare = sum(dr.^2, 2);
end

if isfield(rowCache, 'r_bare') && ~isempty(rowCache.r_bare)
    rBare = double(rowCache.r_bare(:));
else
    rBare = sqrt(r2Bare);
end

if opt.Softening == 0
    if isfield(rowCache, 'inv_r3_bare') && isfield(rowCache, 'inv_r5_bare') && ...
            ~isempty(rowCache.inv_r3_bare) && ~isempty(rowCache.inv_r5_bare)

        invR3 = double(rowCache.inv_r3_bare(:));
        invR5 = double(rowCache.inv_r5_bare(:));
    else
        invR3 = 1 ./ (r2Bare .* rBare);
        invR5 = invR3 ./ r2Bare;
    end
else
    r2 = r2Bare + opt.Softening^2;
    r = sqrt(r2);
    invR3 = 1 ./ (r2 .* r);
    invR5 = invR3 ./ r2;
end

if opt.UseThole
    if isfield(rowCache, 'thole_f3') && isfield(rowCache, 'thole_f5') && ...
            ~isempty(rowCache.thole_f3) && ~isempty(rowCache.thole_f5)

        f3 = double(rowCache.thole_f3(:));
        f5 = double(rowCache.thole_f5(:));
    else
        rowInd = local_row_indices_from_row_ptr(rowPtr);

        activeSites = problem.activeSites(:);
        iFull = activeSites(rowInd);
        jFull = activeSites(colIdx);

        tf = thole.thole_f3f5_factors( ...
            rBare, ...
            sys.site_alpha(iFull), ...
            sys.site_alpha(jFull), ...
            sys.thole_a);

        f3 = tf.f3(:);
        f5 = tf.f5(:);
    end
else
    f3 = ones(nEntries, 1);
    f5 = ones(nEntries, 1);
end

c5 = 3 .* f5 .* invR5;
c3 = f3 .* invR3;

Txx = c5 .* dx .* dx - c3;
Txy = c5 .* dx .* dy;
Txz = c5 .* dx .* dz;

Tyy = c5 .* dy .* dy - c3;
Tyz = c5 .* dy .* dz;

Tzz = c5 .* dz .* dz - c3;

op = struct();
op.mode = 'nonperiodic';
op.kind = 'matrix_free';
op.backend = 'nonperiodic_rowcache_apply';

op.nPolSites = nPol;
op.size = [3*nPol 3*nPol];

op.apply = @(muVec) local_apply_rowcache( ...
    muVec, nPol, rowPtr, colIdx, Txx, Txy, Txz, Tyy, Tyz, Tzz);

op.apply_row = @(rowLocal, muVec) local_apply_rowcache_row( ...
    rowLocal, muVec, nPol, rowPtr, colIdx, Txx, Txy, Txz, Tyy, Tyz, Tzz);

op.row_cache = struct();
op.row_cache.source = 'geom.build_active_row_cache_nonperiodic';
op.row_cache.row_ptr = rowPtr;
op.row_cache.col_idx = colIdx;
op.row_cache.n_entries = nEntries;

% Raw row-cache arrays used by solve_scf_sor fast path.
%
% This mirrors the old matrix-free SOR implementation: the solver performs
% the contraction from dr/f3/f5/invR3/invR5 directly inside the row sweep.
% That is faster for SOR than calling op.apply_row per row or using the
% precomputed tensor-coefficient path.
op.row_cache.dr = dr;
op.row_cache.r2_bare = r2Bare;
op.row_cache.r_bare = rBare;
op.row_cache.inv_r3_bare = invR3;
op.row_cache.inv_r5_bare = invR5;
op.row_cache.thole_f3 = f3;
op.row_cache.thole_f5 = f5;

% Tensor coefficients are retained for op.apply/op.apply_row and tests.
% solve_scf_sor no longer uses these for its production fast path.
op.row_cache.Txx = Txx;
op.row_cache.Txy = Txy;
op.row_cache.Txz = Txz;
op.row_cache.Tyy = Tyy;
op.row_cache.Tyz = Tyz;
op.row_cache.Tzz = Tzz;

op.info = struct();
op.info.assembly_backend = 'nonperiodic_rowcache_apply';
op.info.nPolSites = nPol;
op.info.nEntriesDirected = nEntries;
op.info.rcut = opt.Rcut;
op.info.use_cutoff = true;
op.info.use_mex = opt.UseMex;
op.info.cache_time = cacheTime;
op.info.cache_info = local_summarize_cache(rowCache);

op.capabilities = struct();
op.capabilities.apply = true;
op.capabilities.dense_matrix = false;
op.capabilities.row_update = true;

op.params = struct();
op.params.use_thole = opt.UseThole;
op.params.softening = opt.Softening;
op.params.rcut = opt.Rcut;
op.params.use_mex = opt.UseMex;

end

% =========================================================================
% Apply helpers
% =========================================================================

function Tmu = local_apply_rowcache(muVec, nPol, rowPtr, colIdx, Txx, Txy, Txz, Tyy, Tyz, Tzz)

muVec = muVec(:);

if numel(muVec) ~= 3*nPol
    error('thole:build_nonperiodic_rowcache_operator:BadApplyVectorSize', ...
        'muVec must have length 3*nPolSites.');
end

mux = muVec(1:3:end);
muy = muVec(2:3:end);
muz = muVec(3:3:end);

rowInd = local_row_indices_from_row_ptr(rowPtr);

EentryX = Txx .* mux(colIdx) + Txy .* muy(colIdx) + Txz .* muz(colIdx);
EentryY = Txy .* mux(colIdx) + Tyy .* muy(colIdx) + Tyz .* muz(colIdx);
EentryZ = Txz .* mux(colIdx) + Tyz .* muy(colIdx) + Tzz .* muz(colIdx);

Ex = accumarray(rowInd, EentryX, [nPol 1], @sum, 0);
Ey = accumarray(rowInd, EentryY, [nPol 1], @sum, 0);
Ez = accumarray(rowInd, EentryZ, [nPol 1], @sum, 0);

Tmu = zeros(3*nPol, 1);
Tmu(1:3:end) = Ex;
Tmu(2:3:end) = Ey;
Tmu(3:3:end) = Ez;

end

function TiMu = local_apply_rowcache_row(rowLocal, muVec, nPol, rowPtr, colIdx, Txx, Txy, Txz, Tyy, Tyz, Tzz)

if ~(isnumeric(rowLocal) && isscalar(rowLocal) && rowLocal == round(rowLocal) && ...
        rowLocal >= 1 && rowLocal <= nPol)
    error('thole:build_nonperiodic_rowcache_operator:BadRowIndex', ...
        'rowLocal must be an integer active-site row in 1:nPolSites.');
end

muVec = muVec(:);

if numel(muVec) ~= 3*nPol
    error('thole:build_nonperiodic_rowcache_operator:BadApplyVectorSize', ...
        'muVec must have length 3*nPolSites.');
end

mux = muVec(1:3:end);
muy = muVec(2:3:end);
muz = muVec(3:3:end);

lo = rowPtr(rowLocal);
hi = rowPtr(rowLocal + 1) - 1;

if hi < lo
    TiMu = [0; 0; 0];
    return;
end

idx = lo:hi;
cols = colIdx(idx);

Ex = sum(Txx(idx) .* mux(cols) + Txy(idx) .* muy(cols) + Txz(idx) .* muz(cols));
Ey = sum(Txy(idx) .* mux(cols) + Tyy(idx) .* muy(cols) + Tyz(idx) .* muz(cols));
Ez = sum(Txz(idx) .* mux(cols) + Tyz(idx) .* muy(cols) + Tzz(idx) .* muz(cols));

TiMu = [Ex; Ey; Ez];

end

% =========================================================================
% Cache normalization helpers
% =========================================================================

function rowPtr = local_normalize_row_ptr(rowPtr)

rowPtr = double(rowPtr(:));

if isempty(rowPtr)
    error('thole:build_nonperiodic_rowcache_operator:EmptyRowPtr', ...
        'rowCache.row_ptr is empty.');
end

% Accept either C-style 0-based row_ptr or MATLAB 1-based row_ptr.
if rowPtr(1) == 0
    rowPtr = rowPtr + 1;
end

if rowPtr(1) ~= 1
    error('thole:build_nonperiodic_rowcache_operator:BadRowPtrBase', ...
        'rowCache.row_ptr must start at 0 or 1.');
end

if any(diff(rowPtr) < 0)
    error('thole:build_nonperiodic_rowcache_operator:BadRowPtrMonotonicity', ...
        'rowCache.row_ptr must be nondecreasing.');
end

end

function colIdx = local_normalize_col_idx(colIdx)

colIdx = double(colIdx(:));

if isempty(colIdx)
    return;
end

% Accept either C-style 0-based col_idx or MATLAB 1-based col_idx.
if min(colIdx) == 0
    colIdx = colIdx + 1;
end

end

function colActive = local_map_columns_to_active(colIdxRaw, sys, problem)

nPol = problem.nPolSites;

if isempty(colIdxRaw)
    colActive = colIdxRaw;
    return;
end

% If indices already fit in 1:nPol, treat them as active-local.
if max(colIdxRaw) <= nPol
    colActive = colIdxRaw;

    if any(colActive < 1)
        error('thole:build_nonperiodic_rowcache_operator:BadActiveColIdx', ...
            'Active-local column indices must be positive.');
    end

    return;
end

% Otherwise treat them as full-site indices and map back to active-local.
nSites = size(sys.site_pos, 1);

if any(colIdxRaw < 1) || any(colIdxRaw > nSites)
    error('thole:build_nonperiodic_rowcache_operator:BadFullColIdx', ...
        'Full-site column indices are outside 1:nSites.');
end

fullToActive = zeros(nSites, 1);
fullToActive(problem.activeSites(:)) = 1:nPol;

colActive = fullToActive(colIdxRaw);

if any(colActive == 0)
    error('thole:build_nonperiodic_rowcache_operator:ColOutsideActiveSpace', ...
        'rowCache.col_idx contains a site outside problem.activeSites.');
end

end

function rowInd = local_row_indices_from_row_ptr(rowPtr)

counts = diff(rowPtr(:));
nRows = numel(counts);

rowInd = repelem((1:nRows).', counts);

end

function summary = local_summarize_cache(cache)

summary = struct();

fields = {
    'n_entries'
    'nEntries'
    'nEntriesDirected'
    'rcut'
    'used_mex'
    'use_mex'
    'backend'
    'build_time'
    'timing'
};

for k = 1:numel(fields)
    name = fields{k};

    if isfield(cache, name)
        summary.(name) = cache.(name);
    end
end

if isfield(cache, 'col_idx')
    summary.n_entries_returned = numel(cache.col_idx);
end

end