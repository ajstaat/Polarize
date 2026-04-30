function op = build_periodic_p3m_paircache_operator(sys, problem, opt)
%BUILD_PERIODIC_P3M_PAIRCACHE_OPERATOR Build matrix-free periodic P3M operator.
%
% Private backend builder used by thole.make_polarization_operator.
%
% Returns:
%   op.kind    = 'matrix_free'
%   op.backend = 'periodic_p3m_paircache_apply'
%
% This backend stores:
%   - periodic real-space row cache
%   - periodic P3M reciprocal dipole cache
%   - analytic self/surface blocks
%
% and exposes:
%   op.apply(muVec)
%
% It does not build or store dense op.Tpol and does not expose production
% row-update capability for SOR. Use build_periodic_p3m_rowcache_operator
% for SOR-oriented operators.

io.assert_atomic_units(sys);

if isfield(opt, 'Softening') && opt.Softening ~= 0
    error('thole:build_periodic_p3m_paircache_operator:SofteningUnsupported', ...
        'Periodic P3M matrix-free operators do not support Softening ~= 0.');
end

local_validate_problem(problem);

nPol = problem.nPolSites;
nSites = size(sys.site_pos, 1);

p3mParams = local_get_p3m_params(opt);

rowOpts = struct();
rowOpts.use_mex = local_get_opt(opt, 'UseMex', true);
rowOpts.profile = local_get_opt(opt, 'Profile', false) || local_get_opt(opt, 'Verbose', false);
rowOpts.verbose = local_get_opt(opt, 'Verbose', false);
rowOpts.use_thole = local_get_opt(opt, 'UseThole', true);

ewaldParams = struct();
ewaldParams.alpha = p3mParams.alpha;
ewaldParams.rcut = p3mParams.rcut;
ewaldParams.boundary = p3mParams.boundary;
ewaldParams.kcut = p3mParams.kcut;

tRealCache = tic;
realCacheRaw = geom.build_active_row_cache_periodic(sys, problem, ewaldParams, rowOpts);
realCacheTime = toc(tRealCache);

realCache = local_normalize_real_cache(realCacheRaw, nPol);

tP3MCache = tic;
p3mCache = p3m.build_dipole_cache(sys, local_make_p3m_cache_opts(sys, problem, p3mParams, opt));
p3mCacheTime = toc(tP3MCache);

selfBlock = ewald.self_tensor_block_dipole(p3mParams.alpha);
surfaceBlock = ewald.surface_tensor_block_dipole(p3mCache.lattice, p3mParams.boundary);

op = struct();

op.mode = 'periodic_p3m';
op.kind = 'matrix_free';
op.backend = 'periodic_p3m_paircache_apply';
op.nPolSites = nPol;
op.size = [3*nPol, 3*nPol];

op.apply = @(muVec) local_apply_periodic_p3m( ...
    muVec, nPol, nSites, problem.activeSites(:), ...
    realCache, p3mCache, selfBlock, surfaceBlock);

op.cache = struct();
op.cache.source = 'periodic real-space row cache + periodic P3M dipole cache';
op.cache.nSites = nSites;
op.cache.real_cache = local_summarize_real_cache(realCacheRaw);
op.cache.p3m_cache = local_summarize_p3m_cache(p3mCache);

op.real_cache = realCache;
op.p3m_cache = p3mCache;

op.periodic_p3m_cache = struct();
op.periodic_p3m_cache.real_cache = realCache;
op.periodic_p3m_cache.p3m_cache = p3mCache;
op.periodic_p3m_cache.self_block = selfBlock;
op.periodic_p3m_cache.surface_block = surfaceBlock;
op.periodic_p3m_cache.boundary = p3mParams.boundary;
op.periodic_p3m_cache.alpha = p3mParams.alpha;
op.periodic_p3m_cache.rcut = p3mParams.rcut;
op.periodic_p3m_cache.kcut = p3mParams.kcut;
op.periodic_p3m_cache.mesh_size = p3mParams.mesh_size;
op.periodic_p3m_cache.assignment_order = p3mParams.assignment_order;

% Compatibility convention with periodic Ewald operators.
op.periodic_cache = op.periodic_p3m_cache;

op.info = struct();
op.info.assembly_backend = 'periodic_p3m_paircache_apply';
op.info.nPolSites = nPol;
op.info.nSites = nSites;
op.info.nRealEntriesDirected = realCache.n_entries;
op.info.nK = p3mCache.nK;
op.info.alpha = p3mParams.alpha;
op.info.rcut = p3mParams.rcut;
op.info.kcut = p3mParams.kcut;
op.info.boundary = p3mParams.boundary;
op.info.mesh_size = p3mParams.mesh_size;
op.info.assignment_order = p3mParams.assignment_order;
op.info.derivative_mode = p3mParams.derivative_mode;
op.info.influence_mode = p3mParams.influence_mode;
op.info.fd_stencil = p3mParams.fd_stencil;
op.info.use_cutoff = true;
op.info.use_mex = rowOpts.use_mex;
op.info.real_cache_time = realCacheTime;
op.info.p3m_cache_time = p3mCacheTime;
op.info.cache_time = realCacheTime + p3mCacheTime;
op.info.real_cache_info = local_summarize_real_cache(realCacheRaw);
op.info.p3m_cache_info = local_summarize_p3m_cache(p3mCache);

op.capabilities = struct();
op.capabilities.apply = true;
op.capabilities.dense_matrix = false;
op.capabilities.row_update = false;

op.params = struct();
op.params.use_thole = rowOpts.use_thole;
op.params.softening = 0.0;
op.params.rcut = p3mParams.rcut;
op.params.alpha = p3mParams.alpha;
op.params.kcut = p3mParams.kcut;
op.params.boundary = p3mParams.boundary;
op.params.mesh_size = p3mParams.mesh_size;
op.params.assignment_order = p3mParams.assignment_order;
op.params.derivative_mode = p3mParams.derivative_mode;
op.params.influence_mode = p3mParams.influence_mode;
op.params.fd_stencil = p3mParams.fd_stencil;
op.params.deconvolve_assignment = p3mParams.deconvolve_assignment;
op.params.deconvolution_floor = p3mParams.deconvolution_floor;
op.params.alias_range = p3mParams.alias_range;
op.params.use_mex = rowOpts.use_mex;
end

% =========================================================================
% Apply helper
% =========================================================================

function Tmu = local_apply_periodic_p3m( ...
    muVec, nPol, nSites, activeSites, realCache, p3mCache, selfBlock, surfaceBlock)

muVec = muVec(:);

if numel(muVec) ~= 3*nPol
    error('thole:build_periodic_p3m_paircache_operator:BadApplyVectorSize', ...
        'muVec must have length 3*nPolSites.');
end

mu = zeros(nPol, 3);
mu(:, 1) = muVec(1:3:end);
mu(:, 2) = muVec(2:3:end);
mu(:, 3) = muVec(3:3:end);

E = zeros(nPol, 3);

% Real-space periodic Ewald/Thole contribution.
E = E + thole.apply_periodic_real_cache(realCache, mu);

% Reciprocal-space P3M contribution.
muFull = zeros(nSites, 3);
muFull(activeSites, :) = mu;

[ErecipFull, ~] = p3m.apply_dipole_cache(p3mCache, muFull);
E = E + ErecipFull(activeSites, :);

% Analytic self term.
if ~isempty(selfBlock)
    E = E + mu * selfBlock.';
end

% Surface term.
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

% =========================================================================
% Cache normalization / metadata helpers
% =========================================================================

function realCache = local_normalize_real_cache(cache, nPol)
rowPtr = local_normalize_row_ptr(cache.row_ptr(:));
colIdx = local_normalize_col_idx(cache.col_idx(:));

if numel(rowPtr) ~= nPol + 1
    error('thole:build_periodic_p3m_paircache_operator:BadRowPtr', ...
        'rowCache.row_ptr must have length nPolSites + 1.');
end

nEntries = rowPtr(end) - 1;

if numel(colIdx) ~= nEntries
    error('thole:build_periodic_p3m_paircache_operator:BadColIdxLength', ...
        'numel(rowCache.col_idx) must equal rowPtr(end)-1.');
end

if ~isempty(colIdx) && (min(colIdx) < 1 || max(colIdx) > nPol)
    error('thole:build_periodic_p3m_paircache_operator:BadColIdx', ...
        'Periodic rowCache.col_idx must be active-local indices in 1:nPolSites.');
end

required = {'dr', 'coeff_iso', 'coeff_dyad'};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(cache, name) || isempty(cache.(name))
        error('thole:build_periodic_p3m_paircache_operator:MissingRealCacheField', ...
            'rowCache.%s is required.', name);
    end
end

dr = double(cache.dr);

if size(dr, 1) ~= nEntries || size(dr, 2) ~= 3
    error('thole:build_periodic_p3m_paircache_operator:BadDr', ...
        'rowCache.dr must be nEntries x 3.');
end

coeffIso = double(cache.coeff_iso(:));
coeffDyad = double(cache.coeff_dyad(:));

if numel(coeffIso) ~= nEntries || numel(coeffDyad) ~= nEntries
    error('thole:build_periodic_p3m_paircache_operator:BadCoeffSize', ...
        'rowCache.coeff_iso and coeff_dyad must have length nEntries.');
end

realCache = struct();

realCache.source = 'geom.build_active_row_cache_periodic';
realCache.row_ptr = rowPtr;
realCache.col_idx = colIdx;
realCache.n_entries = nEntries;
realCache.row_ind = local_row_indices_from_row_ptr(rowPtr);

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
    error('thole:build_periodic_p3m_paircache_operator:EmptyRowPtr', ...
        'rowCache.row_ptr is empty.');
end

if rowPtr(1) == 0
    rowPtr = rowPtr + 1;
end

if rowPtr(1) ~= 1
    error('thole:build_periodic_p3m_paircache_operator:BadRowPtrBase', ...
        'rowCache.row_ptr must start at 0 or 1.');
end

if any(diff(rowPtr) < 0)
    error('thole:build_periodic_p3m_paircache_operator:BadRowPtrMonotonicity', ...
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

function summary = local_summarize_p3m_cache(cache)
summary = struct();

fields = {'mode', 'nSites', 'nSources', 'nTargets', 'mesh_size', ...
    'assignment_order', 'alpha', 'derivative_mode', 'influence_mode', ...
    'fd_stencil', 'deconvolve_assignment', 'deconvolution_floor', ...
    'alias_range', 'nK', 'estimated_mesh_gb', 'build_time', ...
    'lattice_convention', 'volume'};

for k = 1:numel(fields)
    name = fields{k};

    if isfield(cache, name)
        summary.(name) = cache.(name);
    end
end
end

function p3mParams = local_get_p3m_params(opt)
p3mParams = struct();

if isfield(opt, 'Ewald') && ~isempty(opt.Ewald)
    ew = opt.Ewald;
else
    ew = struct();
end

p3mParams.alpha = local_get_first_field(opt, ew, {'Alpha', 'alpha'}, []);
p3mParams.kcut = local_get_first_field(opt, ew, {'Kcut', 'kcut'}, NaN);
p3mParams.rcut = local_get_first_field(opt, ew, {'Rcut', 'rcut'}, []);
p3mParams.boundary = local_get_first_field(opt, ew, {'Boundary', 'boundary'}, 'tinfoil');

if isempty(p3mParams.alpha)
    error('thole:build_periodic_p3m_paircache_operator:MissingAlpha', ...
        'Periodic P3M operator requires opt.Alpha or opt.Ewald.alpha.');
end

if isempty(p3mParams.rcut)
    error('thole:build_periodic_p3m_paircache_operator:MissingRcut', ...
        'Periodic P3M operator requires opt.Rcut or opt.Ewald.rcut.');
end

validateattributes(p3mParams.alpha, {'numeric'}, ...
    {'scalar','real','finite','positive'}, mfilename, 'alpha');

validateattributes(p3mParams.rcut, {'numeric'}, ...
    {'scalar','real','finite','positive'}, mfilename, 'rcut');

if ~isnan(p3mParams.kcut)
    validateattributes(p3mParams.kcut, {'numeric'}, ...
        {'scalar','real','finite','positive'}, mfilename, 'kcut');
end

p3mParams.alpha = double(p3mParams.alpha);
p3mParams.kcut = double(p3mParams.kcut);
p3mParams.rcut = double(p3mParams.rcut);

boundary = lower(char(string(p3mParams.boundary)));

if ~ismember(boundary, {'tinfoil','vacuum'})
    error('thole:build_periodic_p3m_paircache_operator:BadBoundary', ...
        'Periodic P3M boundary must be ''tinfoil'' or ''vacuum''.');
end

p3mParams.boundary = boundary;

if ~isfield(opt, 'MeshSize') || isempty(opt.MeshSize)
    error('thole:build_periodic_p3m_paircache_operator:MissingMeshSize', ...
        'Periodic P3M operator requires opt.MeshSize.');
end

meshSize = double(opt.MeshSize(:).');

if numel(meshSize) ~= 3 || any(meshSize <= 0) || any(meshSize ~= round(meshSize))
    error('thole:build_periodic_p3m_paircache_operator:BadMeshSize', ...
        'opt.MeshSize must be a 1x3 positive integer vector.');
end

p3mParams.mesh_size = meshSize;

p3mParams.assignment_order = double(local_get_opt(opt, 'AssignmentOrder', 4));

if p3mParams.assignment_order <= 0 || p3mParams.assignment_order ~= round(p3mParams.assignment_order)
    error('thole:build_periodic_p3m_paircache_operator:BadAssignmentOrder', ...
        'opt.AssignmentOrder must be a positive integer.');
end

p3mParams.derivative_mode = lower(char(string(local_get_opt(opt, 'DerivativeMode', 'spectral'))));
p3mParams.influence_mode = lower(char(string(local_get_opt(opt, 'InfluenceMode', 'ewald'))));
p3mParams.fd_stencil = lower(char(string(local_get_opt(opt, 'FDStencil', 'central2'))));

p3mParams.deconvolve_assignment = logical(local_get_opt(opt, 'DeconvolveAssignment', true));
p3mParams.deconvolution_floor = double(local_get_opt(opt, 'DeconvolutionFloor', 1e-8));
p3mParams.alias_range = double(local_get_opt(opt, 'AliasRange', 2));
end

function p3mOpts = local_make_p3m_cache_opts(sys, problem, p3mParams, opt)
targetMask = false(size(sys.site_pos, 1), 1);
sourceMask = false(size(sys.site_pos, 1), 1);

targetMask(problem.activeSites(:)) = true;
sourceMask(problem.activeSites(:)) = true;

p3mOpts = struct();

p3mOpts.ewald = struct();
p3mOpts.ewald.alpha = p3mParams.alpha;
p3mOpts.ewald.rcut = p3mParams.rcut;
p3mOpts.ewald.boundary = p3mParams.boundary;

if ~isnan(p3mParams.kcut)
    p3mOpts.ewald.kcut = p3mParams.kcut;
end

p3mOpts.mesh_size = p3mParams.mesh_size;
p3mOpts.assignment_order = p3mParams.assignment_order;

p3mOpts.target_mask = targetMask;
p3mOpts.source_mask = sourceMask;

p3mOpts.derivative_mode = p3mParams.derivative_mode;
p3mOpts.influence_mode = p3mParams.influence_mode;
p3mOpts.fd_stencil = p3mParams.fd_stencil;

p3mOpts.deconvolve_assignment = p3mParams.deconvolve_assignment;
p3mOpts.deconvolution_floor = p3mParams.deconvolution_floor;
p3mOpts.alias_range = p3mParams.alias_range;

p3mOpts.verbose = local_get_opt(opt, 'Verbose', false);
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

function local_validate_problem(problem)
required = {'activeSites', 'nPolSites'};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(problem, name) || isempty(problem.(name))
        error('thole:build_periodic_p3m_paircache_operator:MissingProblemField', ...
            'problem.%s is required.', name);
    end
end

if numel(problem.activeSites) ~= problem.nPolSites
    error('thole:build_periodic_p3m_paircache_operator:BadProblemSize', ...
        'numel(problem.activeSites) must equal problem.nPolSites.');
end
end