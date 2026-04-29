function op = build_periodic_rowcache_operator(sys, problem, opt)
%BUILD_PERIODIC_ROWCACHE_OPERATOR Build SOR-oriented periodic Ewald operator.
%
% Private backend builder used by thole.make_polarization_operator.
%
% Returns:
%   op.kind    = 'matrix_free'
%   op.backend = 'periodic_rowcache_apply'
%
% This backend stores:
%   - periodic real-space row cache
%   - periodic reciprocal k-space cache
%   - analytic self/surface blocks
%
% and exposes:
%   op.apply(muVec)
%   op.apply_row(rowLocal, muVec)       slow validation/fallback path only
%
% It sets:
%   op.capabilities.row_update = true
%
% The production SOR solver should not loop over op.apply_row for this
% backend. solve_scf_sor should add a periodic raw fast path that uses
% op.periodic_cache directly, including live/incremental reciprocal source
% sums.

io.assert_atomic_units(sys);

if isfield(opt, 'Softening') && opt.Softening ~= 0
    error('thole:build_periodic_rowcache_operator:SofteningUnsupported', ...
        'Periodic Ewald row-cache operators do not support Softening ~= 0.');
end

nPol = problem.nPolSites;
nSites = size(sys.site_pos, 1);

ewaldParams = local_get_ewald_params(opt);

rowOpts = struct();
rowOpts.use_mex = local_get_opt(opt, 'UseMex', true);
rowOpts.profile = local_get_opt(opt, 'Profile', false) || local_get_opt(opt, 'Verbose', false);
rowOpts.use_thole = local_get_opt(opt, 'UseThole', true);

kOpts = struct();
kOpts.kspace_mode = local_get_opt(opt, 'KspaceMode', 'auto');
kOpts.k_block_size = local_get_opt(opt, 'KBlockSize', 2048);
kOpts.kspace_memory_limit_gb = local_get_opt(opt, 'KspaceMemoryLimitGB', 8);
kOpts.verbose = local_get_opt(opt, 'Verbose', false);

tRealCache = tic;
realCacheRaw = geom.build_active_row_cache_periodic(sys, problem, ewaldParams, rowOpts);
realCacheTime = toc(tRealCache);

realCache = local_normalize_real_cache(realCacheRaw, nPol);

tKCache = tic;
kCache = geom.build_periodic_kspace_cache(sys, problem, ewaldParams, kOpts);
kCacheTime = toc(tKCache);

selfBlock = ewald.self_tensor_block_dipole(ewaldParams.alpha);
surfaceBlock = ewald.surface_tensor_block_dipole(kCache.lattice, ewaldParams.boundary);

op = struct();

op.mode = 'periodic_ewald';
op.kind = 'matrix_free';
op.backend = 'periodic_rowcache_apply';

op.nPolSites = nPol;
op.size = [3*nPol, 3*nPol];

op.apply = @(muVec) local_apply_periodic_ewald( ...
    muVec, nPol, realCache, kCache, selfBlock, surfaceBlock);

% Slow correctness/fallback row apply. Production SOR should bypass this.
op.apply_row = @(rowLocal, muVec) local_apply_periodic_ewald_row_slow( ...
    rowLocal, muVec, nPol, realCache, kCache, selfBlock, surfaceBlock);

op.row_cache = realCache;
op.k_cache = kCache;

op.periodic_cache = struct();
op.periodic_cache.real_cache = realCache;
op.periodic_cache.k_cache = kCache;
op.periodic_cache.self_block = selfBlock;
op.periodic_cache.surface_block = surfaceBlock;
op.periodic_cache.self_coeff = selfBlock(1, 1);
op.periodic_cache.surface_coeff = surfaceBlock(1, 1);
op.periodic_cache.boundary = ewaldParams.boundary;
op.periodic_cache.alpha = ewaldParams.alpha;
op.periodic_cache.rcut = ewaldParams.rcut;
op.periodic_cache.kcut = ewaldParams.kcut;
op.periodic_cache.use_mex_kspace = local_get_opt(opt, 'UseMexKspace', true);

op.cache = struct();
op.cache.source = 'periodic real-space row cache + periodic k-space cache';
op.cache.nSites = nSites;
op.cache.real_cache = local_summarize_real_cache(realCacheRaw);
op.cache.k_cache = local_summarize_k_cache(kCache);

op.info = struct();
op.info.assembly_backend = 'periodic_rowcache_apply';
op.info.nPolSites = nPol;
op.info.nSites = nSites;
op.info.nRealEntriesDirected = realCache.n_entries;
op.info.nK = kCache.num_kvec;
op.info.alpha = ewaldParams.alpha;
op.info.rcut = ewaldParams.rcut;
op.info.kcut = ewaldParams.kcut;
op.info.boundary = ewaldParams.boundary;
op.info.use_cutoff = true;
op.info.use_mex = rowOpts.use_mex;
op.info.real_cache_time = realCacheTime;
op.info.k_cache_time = kCacheTime;
op.info.cache_time = realCacheTime + kCacheTime;
op.info.real_cache_info = local_summarize_real_cache(realCacheRaw);
op.info.k_cache_info = local_summarize_k_cache(kCache);

op.capabilities = struct();
op.capabilities.apply = true;
op.capabilities.dense_matrix = false;
op.capabilities.row_update = true;

op.params = struct();
op.params.use_thole = rowOpts.use_thole;
op.params.softening = 0.0;
op.params.rcut = ewaldParams.rcut;
op.params.alpha = ewaldParams.alpha;
op.params.kcut = ewaldParams.kcut;
op.params.boundary = ewaldParams.boundary;
op.params.use_mex = rowOpts.use_mex;
op.params.kspace_mode = kOpts.kspace_mode;
op.params.k_block_size = kOpts.k_block_size;
op.params.kspace_memory_limit_gb = kOpts.kspace_memory_limit_gb;
end

% =========================================================================
% Apply helpers
% =========================================================================

function TiMu = local_apply_periodic_ewald_row_slow(rowLocal, muVec, nPol, realCache, kCache, selfBlock, surfaceBlock)
if ~(isnumeric(rowLocal) && isscalar(rowLocal) && rowLocal == round(rowLocal) && ...
        rowLocal >= 1 && rowLocal <= nPol)
    error('thole:build_periodic_rowcache_operator:BadRowIndex', ...
        'rowLocal must be an integer active-site row in 1:nPolSites.');
end

Tmu = local_apply_periodic_ewald(muVec, nPol, realCache, kCache, selfBlock, surfaceBlock);

idx = 3*(rowLocal - 1) + (1:3);
TiMu = Tmu(idx);
end

function Tmu = local_apply_periodic_ewald(muVec, nPol, realCache, kCache, selfBlock, surfaceBlock)
muVec = muVec(:);

if numel(muVec) ~= 3*nPol
    error('thole:build_periodic_rowcache_operator:BadApplyVectorSize', ...
        'muVec must have length 3*nPolSites.');
end

mu = zeros(nPol, 3);
mu(:, 1) = muVec(1:3:end);
mu(:, 2) = muVec(2:3:end);
mu(:, 3) = muVec(3:3:end);

E = zeros(nPol, 3);

% Real-space contribution.
if realCache.n_entries > 0
    rowInd = local_row_indices_from_row_ptr(realCache.row_ptr);
    cols = realCache.col_idx;

    muSrc = mu(cols, :);
    dr = realCache.dr;

    muDotR = sum(muSrc .* dr, 2);

    Eentry = realCache.coeff_iso .* muSrc + ...
             realCache.coeff_dyad .* muDotR .* dr;

    E(:, 1) = E(:, 1) + accumarray(rowInd, Eentry(:, 1), [nPol 1], @sum, 0);
    E(:, 2) = E(:, 2) + accumarray(rowInd, Eentry(:, 2), [nPol 1], @sum, 0);
    E(:, 3) = E(:, 3) + accumarray(rowInd, Eentry(:, 3), [nPol 1], @sum, 0);
end

% Reciprocal-space contribution.
E = E + local_apply_kspace(mu, kCache);

% Self.
if ~isempty(selfBlock)
    E = E + mu * selfBlock.';
end

% Surface.
if ~isempty(surfaceBlock) && any(surfaceBlock(:) ~= 0)
    M = sum(mu, 1);
    Esurf = (surfaceBlock * M.').';
    E = E + repmat(Esurf, nPol, 1);
end

Tmu = zeros(3*nPol, 1);
Tmu(1:3:end) = E(:, 1);
Tmu(2:3:end) = E(:, 2);
Tmu(3:3:end) = E(:, 3);
end

function Erecip = local_apply_kspace(mu, kCache)
nPol = size(mu, 1);
Erecip = zeros(nPol, 3);

if ~isfield(kCache, 'num_kvec') || kCache.num_kvec == 0
    return;
end

storageMode = 'full';
if isfield(kCache, 'storage_mode') && ~isempty(kCache.storage_mode)
    storageMode = lower(char(string(kCache.storage_mode)));
end

switch storageMode
    case 'full'
        kvecs = kCache.kvecs;
        kvecsT = kCache.kvecs_T;
        twoPref = kCache.two_pref(:);

        if isfield(kCache, 'cos_phase') && ~isempty(kCache.cos_phase)
            cosPhase = kCache.cos_phase;
            sinPhase = kCache.sin_phase;
        else
            phase = kCache.active_pos * kvecsT;
            cosPhase = cos(phase);
            sinPhase = sin(phase);
        end

        v = mu * kvecsT;

        A = sum(cosPhase .* v, 1);
        B = sum(sinPhase .* v, 1);

        phaseFactor = cosPhase .* A + sinPhase .* B;
        W = phaseFactor .* twoPref.';

        Erecip = Erecip + W * kvecs;

    case 'blocked'
        if ~isfield(kCache, 'blocks') || isempty(kCache.blocks)
            error('thole:build_periodic_rowcache_operator:MissingKBlocks', ...
                'Blocked k-space cache has no blocks.');
        end

        for b = 1:numel(kCache.blocks)
            blk = kCache.blocks(b);

            if blk.nk == 0
                continue;
            end

            kb = blk.kvecs;
            kbT = blk.kvecs_T;
            twoPref = blk.two_pref(:);

            if isfield(blk, 'cos_phase') && ~isempty(blk.cos_phase) && ...
                    size(blk.cos_phase, 2) == blk.nk
                cosPhase = blk.cos_phase;
                sinPhase = blk.sin_phase;
            else
                phase = kCache.active_pos * kbT;
                cosPhase = cos(phase);
                sinPhase = sin(phase);
            end

            v = mu * kbT;

            A = sum(cosPhase .* v, 1);
            B = sum(sinPhase .* v, 1);

            phaseFactor = cosPhase .* A + sinPhase .* B;
            W = phaseFactor .* twoPref.';

            Erecip = Erecip + W * kb;
        end

    otherwise
        error('thole:build_periodic_rowcache_operator:BadKspaceStorageMode', ...
            'Unsupported kCache.storage_mode "%s".', storageMode);
end
end

% =========================================================================
% Cache normalization / metadata helpers
% =========================================================================

function realCache = local_normalize_real_cache(cache, nPol)
rowPtr = local_normalize_row_ptr(cache.row_ptr(:));
colIdx = local_normalize_col_idx(cache.col_idx(:));

if numel(rowPtr) ~= nPol + 1
    error('thole:build_periodic_rowcache_operator:BadRowPtr', ...
        'rowCache.row_ptr must have length nPolSites + 1.');
end

nEntries = rowPtr(end) - 1;

if numel(colIdx) ~= nEntries
    error('thole:build_periodic_rowcache_operator:BadColIdxLength', ...
        'numel(rowCache.col_idx) must equal rowPtr(end)-1.');
end

if ~isempty(colIdx) && (min(colIdx) < 1 || max(colIdx) > nPol)
    error('thole:build_periodic_rowcache_operator:BadColIdx', ...
        'Periodic rowCache.col_idx must be active-local indices in 1:nPolSites.');
end

required = {'dr', 'coeff_iso', 'coeff_dyad'};
for k = 1:numel(required)
    name = required{k};
    if ~isfield(cache, name) || isempty(cache.(name))
        error('thole:build_periodic_rowcache_operator:MissingRealCacheField', ...
            'rowCache.%s is required.', name);
    end
end

dr = double(cache.dr);
if size(dr, 1) ~= nEntries || size(dr, 2) ~= 3
    error('thole:build_periodic_rowcache_operator:BadDr', ...
        'rowCache.dr must be nEntries x 3.');
end

coeffIso = double(cache.coeff_iso(:));
coeffDyad = double(cache.coeff_dyad(:));

if numel(coeffIso) ~= nEntries || numel(coeffDyad) ~= nEntries
    error('thole:build_periodic_rowcache_operator:BadCoeffSize', ...
        'rowCache.coeff_iso and coeff_dyad must have length nEntries.');
end

realCache = struct();
realCache.source = 'geom.build_active_row_cache_periodic';
realCache.row_ptr = rowPtr;
realCache.col_idx = colIdx;
realCache.n_entries = nEntries;

realCache.dr = dr;
realCache.coeff_iso = coeffIso;
realCache.coeff_dyad = coeffDyad;

if isfield(cache, 'r2_bare')
    realCache.r2_bare = double(cache.r2_bare(:));
else
    realCache.r2_bare = sum(dr.^2, 2);
end

if isfield(cache, 'r_bare')
    realCache.r_bare = double(cache.r_bare(:));
else
    realCache.r_bare = sqrt(realCache.r2_bare);
end

copyFields = {'activeSites', 'nActive', 'alpha', 'rcut', 'use_thole', ...
    'nPairsUndirected', 'nEntriesDirected', 'nInteractions', ...
    'isPeriodic', 'lattice', 'H', 'G', 'V', 'Lmin', ...
    'lattice_convention', 'convention', 'backend', 'profileTimes'};

for k = 1:numel(copyFields)
    name = copyFields{k};
    if isfield(cache, name)
        realCache.(name) = cache.(name);
    end
end
end

function rowPtr = local_normalize_row_ptr(rowPtr)
rowPtr = double(rowPtr(:));

if isempty(rowPtr)
    error('thole:build_periodic_rowcache_operator:EmptyRowPtr', ...
        'rowCache.row_ptr is empty.');
end

if rowPtr(1) == 0
    rowPtr = rowPtr + 1;
end

if rowPtr(1) ~= 1
    error('thole:build_periodic_rowcache_operator:BadRowPtrBase', ...
        'rowCache.row_ptr must start at 0 or 1.');
end

if any(diff(rowPtr) < 0)
    error('thole:build_periodic_rowcache_operator:BadRowPtrMonotonicity', ...
        'rowCache.row_ptr must be nondecreasing.');
end
end

function colIdx = local_normalize_col_idx(colIdx)
colIdx = double(colIdx(:));

if isempty(colIdx)
    return;
end

if min(colIdx) == 0
    colIdx = colIdx + 1;
end
end

function rowInd = local_row_indices_from_row_ptr(rowPtr)
counts = diff(rowPtr(:));
nRows = numel(counts);
rowInd = repelem((1:nRows).', counts);
end

function summary = local_summarize_real_cache(cache)
summary = struct();

fields = {'mode', 'nActive', 'nPairsUndirected', 'nEntriesDirected', ...
    'nInteractions', 'alpha', 'rcut', 'use_thole', 'backend', ...
    'Lmin', 'lattice_convention', 'convention'};

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

function summary = local_summarize_k_cache(kCache)
summary = struct();

fields = {'mode', 'num_kvec', 'hkmax', 'alpha', 'kcut', 'boundary', ...
    'storage_mode', 'phase_storage', 'estimated_full_gb', ...
    'memory_limit_gb', 'k_block_size', 'num_blocks', ...
    'convention', 'lattice_convention'};

for k = 1:numel(fields)
    name = fields{k};
    if isfield(kCache, name)
        summary.(name) = kCache.(name);
    end
end

if isfield(kCache, 'kmeta')
    summary.kmeta = kCache.kmeta;
end
end

function ewaldParams = local_get_ewald_params(opt)
ewaldParams = struct();

if isfield(opt, 'Ewald') && ~isempty(opt.Ewald)
    ew = opt.Ewald;
else
    ew = struct();
end

ewaldParams.alpha = local_get_first_field(opt, ew, {'Alpha', 'alpha'}, []);
ewaldParams.kcut = local_get_first_field(opt, ew, {'Kcut', 'kcut'}, []);
ewaldParams.rcut = local_get_first_field(opt, ew, {'Rcut', 'rcut'}, []);

if isempty(ewaldParams.alpha)
    error('thole:build_periodic_rowcache_operator:MissingAlpha', ...
        'Periodic Ewald operator requires opt.Alpha or opt.Ewald.alpha.');
end

if isempty(ewaldParams.kcut)
    error('thole:build_periodic_rowcache_operator:MissingKcut', ...
        'Periodic Ewald operator requires opt.Kcut or opt.Ewald.kcut.');
end

if isempty(ewaldParams.rcut)
    error('thole:build_periodic_rowcache_operator:MissingRcut', ...
        'Periodic Ewald operator requires opt.Rcut or opt.Ewald.rcut.');
end

validateattributes(ewaldParams.alpha, {'numeric'}, ...
    {'scalar','real','finite','positive'}, mfilename, 'alpha');
validateattributes(ewaldParams.kcut, {'numeric'}, ...
    {'scalar','real','finite','positive'}, mfilename, 'kcut');
validateattributes(ewaldParams.rcut, {'numeric'}, ...
    {'scalar','real','finite','positive'}, mfilename, 'rcut');

ewaldParams.alpha = double(ewaldParams.alpha);
ewaldParams.kcut = double(ewaldParams.kcut);
ewaldParams.rcut = double(ewaldParams.rcut);

boundary = local_get_first_field(opt, ew, {'Boundary', 'boundary'}, 'tinfoil');
boundary = lower(char(string(boundary)));

if ~ismember(boundary, {'tinfoil','vacuum'})
    error('thole:build_periodic_rowcache_operator:BadBoundary', ...
        'Periodic Ewald boundary must be ''tinfoil'' or ''vacuum''.');
end

ewaldParams.boundary = boundary;
end

function value = local_get_first_field(opt, ew, names, defaultValue)
value = defaultValue;

for k = 1:numel(names)
    name = names{k};
    if isfield(opt, name) && ~isempty(opt.(name))
        value = opt.(name);
        return;
    end
end

for k = 1:numel(names)
    name = names{k};
    if isfield(ew, name) && ~isempty(ew.(name))
        value = ew.(name);
        return;
    end
end
end

function value = local_get_opt(opt, name, defaultValue)
if isfield(opt, name) && ~isempty(opt.(name))
    value = opt.(name);
else
    value = defaultValue;
end
end