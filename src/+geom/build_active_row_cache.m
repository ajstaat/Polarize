function rowCache = build_active_row_cache(sys, problem, opts, spatial)
%BUILD_ACTIVE_ROW_CACHE Build directed active-space row cache.
%
% Uses MEX accelerator when available for finite-rcut nonperiodic cell-list
% builds, otherwise falls back to the stable MATLAB implementation.
%
% rowCache = geom.build_active_row_cache(sys, problem)
% rowCache = geom.build_active_row_cache(sys, problem, opts)
% rowCache = geom.build_active_row_cache(sys, problem, opts, spatial)
%
% Inputs
%   sys      canonical polarization system
%   problem  struct with field .activeSites
%   opts     optional struct with fields:
%              .rcut
%              .profile
%              .use_mex
%   spatial  optional prebuilt spatial index on active-site positions
%
% Output
%   rowCache struct with fields used by solve_scf_iterative_sor

narginchk(2, 4);
if nargin < 3 || isempty(opts)
    opts = struct();
end

% Default: use MEX automatically if available.
useMex = true;
if isfield(opts, 'use_mex') && ~isempty(opts.use_mex)
    useMex = logical(opts.use_mex);
end

doProfile = false;
if isfield(opts, 'profile') && ~isempty(opts.profile)
    doProfile = logical(opts.profile);
end

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
    error('geom:build_active_row_cache:MissingSitePos', ...
        'sys.site_pos is required and may not be empty.');
end
if ~isfield(problem, 'activeSites') || isempty(problem.activeSites)
    error('geom:build_active_row_cache:MissingActiveSites', ...
        'problem.activeSites is required and may not be empty.');
end

pos = sys.site_pos;
nSites = size(pos, 1);

activeSites = problem.activeSites(:);
if any(activeSites > nSites)
    error('geom:build_active_row_cache:ActiveSiteOutOfRange', ...
        'problem.activeSites contains indices outside sys.site_pos.');
end
if numel(unique(activeSites)) ~= numel(activeSites)
    error('geom:build_active_row_cache:DuplicateActiveSites', ...
        'problem.activeSites must not contain duplicate indices.');
end

nActive = numel(activeSites);
posAct = pos(activeSites, :);

rcut = inf;
if isfield(opts, 'rcut') && ~isempty(opts.rcut)
    rcut = opts.rcut;
end
if ~(isscalar(rcut) && isreal(rcut) && isfinite(rcut) && rcut > 0) && ~isinf(rcut)
    error('geom:build_active_row_cache:BadRcut', ...
        'opts.rcut must be a positive scalar or Inf.');
end

profileTimes = struct();
tTotal = tic;

% -------------------------------------------------------------------------
% Spatial index
% -------------------------------------------------------------------------
if nargin < 4 || isempty(spatial)
    spatialOpts = struct();
    spatialOpts.isPeriodic = false;

    if isfinite(rcut)
        spatialOpts.method = 'auto';
        spatialOpts.cutoff = rcut;
    else
        spatialOpts.method = 'bruteforce';
    end

    t0 = tic;
    spatial = geom.build_spatial_index(posAct, spatialOpts);
    profileTimes.build_spatial_index = toc(t0);
else
    profileTimes.build_spatial_index = 0.0;
end

mexAvailable = (exist(['mex_build_active_row_cache_nonperiodic.' mexext], 'file') == 3) || ...
               (exist('mex_build_active_row_cache_nonperiodic', 'file') == 3);

useMexPath = useMex && mexAvailable && ...
    isfinite(rcut) && ...
    isstruct(spatial) && ...
    isfield(spatial, 'backend') && strcmp(spatial.backend, 'cell_list') && ...
    isfield(spatial, 'isPeriodic') && ~spatial.isPeriodic;

if useMexPath
    % ---------------------------------------------------------------------
    % Prepare optional Thole inputs on active sites
    % ---------------------------------------------------------------------
    alphaAct = [];
    tholeAAct = [];

    if isfield(sys, 'site_alpha') && ~isempty(sys.site_alpha) && ...
       isfield(sys, 'thole_a') && ~isempty(sys.thole_a)

        alphaFull = sys.site_alpha(:);
        if numel(alphaFull) ~= nSites
            error('geom:build_active_row_cache:BadSiteAlpha', ...
                'sys.site_alpha must have length equal to size(sys.site_pos,1).');
        end
        alphaAct = alphaFull(activeSites);

        a = sys.thole_a;
        if isscalar(a)
            tholeAAct = double(a);
        else
            a = a(:);
            if numel(a) ~= nSites
                error('geom:build_active_row_cache:BadTholeA', ...
                    'Non-scalar sys.thole_a must have length equal to size(sys.site_pos,1).');
            end
            tholeAAct = a(activeSites);
        end
    end

    % ---------------------------------------------------------------------
    % MEX build
    % ---------------------------------------------------------------------
    t0 = tic;
    [row_ptr, col_idx, dr, r2_bare, r_bare, inv_r3_bare, inv_r5_bare, ...
        thole_f3, thole_f5, nPairsUndirected] = ...
        mex_build_active_row_cache_nonperiodic( ...
            spatial.pos, ...
            double(spatial.grid_shape), ...
            double(spatial.bin_head), ...
            double(spatial.bin_next), ...
            double(spatial.neighbor_offsets), ...
            double(rcut), ...
            alphaAct, ...
            tholeAAct);
    profileTimes.query_pairs = toc(t0);

    profileTimes.undirected_geometry = 0.0;
    profileTimes.thole_undirected = 0.0;
    profileTimes.csr_counts = 0.0;
    profileTimes.csr_fill = 0.0;
    profileTimes.total = toc(tTotal);

    rowCache = struct();
    rowCache.activeSites = activeSites;
    rowCache.nActive = nActive;
    rowCache.row_ptr = row_ptr;
    rowCache.col_idx = col_idx;
    rowCache.dr = dr;
    rowCache.r2_bare = r2_bare;
    rowCache.r_bare = r_bare;
    rowCache.inv_r3_bare = inv_r3_bare;
    rowCache.inv_r5_bare = inv_r5_bare;
    rowCache.rcut = rcut;
    rowCache.nPairsUndirected = nPairsUndirected;
    rowCache.nEntriesDirected = numel(col_idx);

    if ~isempty(thole_f3)
        rowCache.thole_f3 = thole_f3;
        rowCache.thole_f5 = thole_f5;
    end

    if doProfile
        rowCache.profileTimes = profileTimes;
        fprintf('build_active_row_cache timing summary:\n');
        fprintf('  build_spatial_index   : %.6f s\n', profileTimes.build_spatial_index);
        fprintf('  query_pairs           : %.6f s\n', profileTimes.query_pairs);
        fprintf('  undirected_geometry   : %.6f s\n', profileTimes.undirected_geometry);
        fprintf('  thole_undirected      : %.6f s\n', profileTimes.thole_undirected);
        fprintf('  csr_counts            : %.6f s\n', profileTimes.csr_counts);
        fprintf('  csr_fill              : %.6f s\n', profileTimes.csr_fill);
        fprintf('  total                 : %.6f s\n', profileTimes.total);
        fprintf('  nPairsUndirected      : %d\n', rowCache.nPairsUndirected);
        fprintf('  nEntriesDirected      : %d\n', rowCache.nEntriesDirected);
        fprintf('  backend               : MEX\n');
    end
    return;
end

% -------------------------------------------------------------------------
% MATLAB fallback
% -------------------------------------------------------------------------
rowCache = local_build_active_row_cache_matlab_stable(sys, problem, opts, spatial);
end

% =========================================================================
% Stable MATLAB fallback
% =========================================================================
function rowCache = local_build_active_row_cache_matlab_stable(sys, problem, opts, spatial)

doProfile = false;
if isfield(opts, 'profile') && ~isempty(opts.profile)
    doProfile = logical(opts.profile);
end

tTotal = tic;

pos = sys.site_pos;
nSites = size(pos, 1);
activeSites = problem.activeSites(:);
nActive = numel(activeSites);
rcut = inf;
if isfield(opts, 'rcut') && ~isempty(opts.rcut)
    rcut = opts.rcut;
end

posAct = pos(activeSites, :);

profileTimes = struct();

if nargin < 4 || isempty(spatial) || isempty(spatial)
    spatialOpts = struct();
    spatialOpts.isPeriodic = false;
    if isfinite(rcut)
        spatialOpts.method = 'auto';
        spatialOpts.cutoff = rcut;
    else
        spatialOpts.method = 'bruteforce';
    end

    t0 = tic;
    spatial = geom.build_spatial_index(posAct, spatialOpts);
    profileTimes.build_spatial_index = toc(t0);
else
    profileTimes.build_spatial_index = 0.0;
end

queryOpts = struct();
queryOpts.return_r = true;
queryOpts.return_dr = true;
queryOpts.return_full_idx = false;

t0 = tic;
pairs = geom.query_pairs_within_cutoff(spatial, rcut, queryOpts);
profileTimes.query_pairs = toc(t0);

t0 = tic;
pairI = pairs.i(:);
pairJ = pairs.j(:);
drUndir = pairs.dr;
r2Undir = pairs.r2(:);
rUndir = pairs.r(:);

if any(r2Undir <= 0)
    error('geom:build_active_row_cache:NonPositiveDistance', ...
        'Encountered non-positive pair distance while building active row cache.');
end

nPairs = numel(pairI);
nDir = 2 * nPairs;

invR3Undir = zeros(nPairs, 1);
invR5Undir = zeros(nPairs, 1);
if nPairs > 0
    invRUndir = 1 ./ rUndir;
    invR3Undir = invRUndir ./ r2Undir;
    invR5Undir = invR3Undir ./ r2Undir;
end
profileTimes.undirected_geometry = toc(t0);

t0 = tic;
haveThole = false;
f3Undir = [];
f5Undir = [];

if isfield(sys, 'site_alpha') && ~isempty(sys.site_alpha) && ...
   isfield(sys, 'thole_a') && ~isempty(sys.thole_a)

    alpha = sys.site_alpha(:);
    if numel(alpha) ~= nSites
        error('geom:build_active_row_cache:BadSiteAlpha', ...
            'sys.site_alpha must have length equal to size(sys.site_pos,1).');
    end

    haveThole = true;
    if nPairs > 0
        alpha_i = alpha(activeSites(pairI));
        alpha_j = alpha(activeSites(pairJ));
        tf = thole.thole_f3f5_factors(rUndir, alpha_i, alpha_j, sys.thole_a);
        f3Undir = tf.f3;
        f5Undir = tf.f5;
    else
        f3Undir = zeros(0, 1);
        f5Undir = zeros(0, 1);
    end
end
profileTimes.thole_undirected = toc(t0);

t0 = tic;
rowCounts = accumarray([pairI; pairJ], 1, [nActive, 1], @sum, 0);

row_ptr = zeros(nActive + 1, 1);
row_ptr(1) = 1;
if nActive > 0
    row_ptr(2:end) = 1 + cumsum(rowCounts);
end
profileTimes.csr_counts = toc(t0);

t0 = tic;
col_idx = zeros(nDir, 1);
drDir = zeros(nDir, 3);
r2Dir = zeros(nDir, 1);
rDir = zeros(nDir, 1);
invR3Dir = zeros(nDir, 1);
invR5Dir = zeros(nDir, 1);

if haveThole
    f3Dir = zeros(nDir, 1);
    f5Dir = zeros(nDir, 1);
else
    f3Dir = [];
    f5Dir = [];
end

nextPtr = row_ptr(1:end-1);

for p = 1:nPairs
    i = pairI(p);
    j = pairJ(p);

    k = nextPtr(i);
    col_idx(k) = j;
    drDir(k, :) = drUndir(p, :);
    r2Dir(k) = r2Undir(p);
    rDir(k) = rUndir(p);
    invR3Dir(k) = invR3Undir(p);
    invR5Dir(k) = invR5Undir(p);
    if haveThole
        f3Dir(k) = f3Undir(p);
        f5Dir(k) = f5Undir(p);
    end
    nextPtr(i) = k + 1;

    k = nextPtr(j);
    col_idx(k) = i;
    drDir(k, :) = -drUndir(p, :);
    r2Dir(k) = r2Undir(p);
    rDir(k) = rUndir(p);
    invR3Dir(k) = invR3Undir(p);
    invR5Dir(k) = invR5Undir(p);
    if haveThole
        f3Dir(k) = f3Undir(p);
        f5Dir(k) = f5Undir(p);
    end
    nextPtr(j) = k + 1;
end
profileTimes.csr_fill = toc(t0);

rowCache = struct();
rowCache.activeSites = activeSites;
rowCache.nActive = nActive;
rowCache.row_ptr = row_ptr;
rowCache.col_idx = col_idx;
rowCache.dr = drDir;
rowCache.r2_bare = r2Dir;
rowCache.r_bare = rDir;
rowCache.inv_r3_bare = invR3Dir;
rowCache.inv_r5_bare = invR5Dir;
rowCache.rcut = rcut;
rowCache.nPairsUndirected = nPairs;
rowCache.nEntriesDirected = nDir;

if haveThole
    rowCache.thole_f3 = f3Dir;
    rowCache.thole_f5 = f5Dir;
end

profileTimes.total = toc(tTotal);

if doProfile
    rowCache.profileTimes = profileTimes;
    fprintf('build_active_row_cache timing summary:\n');
    fprintf('  build_spatial_index   : %.6f s\n', profileTimes.build_spatial_index);
    fprintf('  query_pairs           : %.6f s\n', profileTimes.query_pairs);
    fprintf('  undirected_geometry   : %.6f s\n', profileTimes.undirected_geometry);
    fprintf('  thole_undirected      : %.6f s\n', profileTimes.thole_undirected);
    fprintf('  csr_counts            : %.6f s\n', profileTimes.csr_counts);
    fprintf('  csr_fill              : %.6f s\n', profileTimes.csr_fill);
    fprintf('  total                 : %.6f s\n', profileTimes.total);
    fprintf('  nPairsUndirected      : %d\n', nPairs);
    fprintf('  nEntriesDirected      : %d\n', nDir);
    fprintf('  backend               : MATLAB\n');
end
end