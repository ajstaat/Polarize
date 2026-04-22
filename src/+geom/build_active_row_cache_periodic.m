function rowCache = build_active_row_cache_periodic(sys, problem, ewaldParams, opts, spatial)
%BUILD_ACTIVE_ROW_CACHE_PERIODIC Build directed active-space row cache
% for periodic real-space Ewald interactions in the single-image regime.
%
% Uses MEX accelerator when available for finite-rcut periodic cell-list
% builds, otherwise falls back to the stable MATLAB implementation.
%
% rowCache = geom.build_active_row_cache_periodic(sys, problem, ewaldParams)
% rowCache = geom.build_active_row_cache_periodic(sys, problem, ewaldParams, opts)
% rowCache = geom.build_active_row_cache_periodic(sys, problem, ewaldParams, opts, spatial)

narginchk(3, 5);
if nargin < 4 || isempty(opts)
    opts = struct();
end

useMex = true;
if isfield(opts, 'use_mex') && ~isempty(opts.use_mex)
    useMex = logical(opts.use_mex);
end

doProfile = false;
if isfield(opts, 'profile') && ~isempty(opts.profile)
    doProfile = logical(opts.profile);
end

useThole = true;
if isfield(opts, 'use_thole') && ~isempty(opts.use_thole)
    useThole = logical(opts.use_thole);
end

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
    error('geom:build_active_row_cache_periodic:MissingSitePos', ...
        'sys.site_pos is required and may not be empty.');
end
if ~isfield(problem, 'activeSites') || isempty(problem.activeSites)
    error('geom:build_active_row_cache_periodic:MissingActiveSites', ...
        'problem.activeSites is required and may not be empty.');
end
if ~isfield(ewaldParams, 'alpha') || isempty(ewaldParams.alpha)
    error('geom:build_active_row_cache_periodic:MissingAlpha', ...
        'ewaldParams.alpha is required.');
end
if ~isfield(ewaldParams, 'rcut') || isempty(ewaldParams.rcut)
    error('geom:build_active_row_cache_periodic:MissingRcut', ...
        'ewaldParams.rcut is required.');
end

alphaEwald = ewaldParams.alpha;
rcut = ewaldParams.rcut;

validateattributes(alphaEwald, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'ewaldParams.alpha');
validateattributes(rcut, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'ewaldParams.rcut');

pos = sys.site_pos;
nSites = size(pos, 1);

activeSites = problem.activeSites(:);
if any(activeSites > nSites)
    error('geom:build_active_row_cache_periodic:ActiveSiteOutOfRange', ...
        'problem.activeSites contains indices outside sys.site_pos.');
end
if numel(unique(activeSites)) ~= numel(activeSites)
    error('geom:build_active_row_cache_periodic:DuplicateActiveSites', ...
        'problem.activeSites must not contain duplicate indices.');
end

nActive = numel(activeSites);
posAct = pos(activeSites, :);

H = local_get_periodic_lattice(sys);
Lmin = geom.shortest_lattice_translation(H);

tol = 1e-12 * max(1, Lmin);
rcutMaxSafe = 0.5 * Lmin - tol;
if ~(rcut < rcutMaxSafe)
    error('geom:build_active_row_cache_periodic:RcutTooLarge', ...
        ['Periodic real-space row cache assumes single-image treatment.\n' ...
         'Require rcut < Lmin/2.\n' ...
         '  Lmin          = %.16g\n' ...
         '  Lmin/2        = %.16g\n' ...
         '  max safe rcut = %.16g\n' ...
         '  current rcut  = %.16g'], ...
         Lmin, 0.5 * Lmin, rcutMaxSafe, rcut);
end

profileTimes = struct();
tTotal = tic;

if nargin < 5 || isempty(spatial)
    spatialOpts = struct();
    spatialOpts.isPeriodic = true;
    spatialOpts.cell = H;
    spatialOpts.method = 'auto';
    spatialOpts.cutoff = rcut;

    t0 = tic;
    spatial = geom.build_spatial_index(posAct, spatialOpts);
    profileTimes.build_spatial_index = toc(t0);
else
    profileTimes.build_spatial_index = 0.0;
end

if ~isstruct(spatial) || ~isfield(spatial, 'backend') || ...
   ~isfield(spatial, 'isPeriodic') || ~spatial.isPeriodic
    error('geom:build_active_row_cache_periodic:BadSpatialIndex', ...
        'Periodic row cache requires a periodic spatial index.');
end

mexAvailable = (exist(['mex_build_active_row_cache_periodic.' mexext], 'file') == 3) || ...
               (exist('mex_build_active_row_cache_periodic', 'file') == 3);

useMexPath = useMex && mexAvailable && ...
    isfinite(rcut) && ...
    isstruct(spatial) && ...
    isfield(spatial, 'backend') && strcmp(spatial.backend, 'cell_list') && ...
    isfield(spatial, 'isPeriodic') && spatial.isPeriodic;

if useMexPath
    alphaAct = [];
    tholeAAct = [];

    if useThole
        if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
            error('geom:build_active_row_cache_periodic:MissingSiteAlpha', ...
                'sys.site_alpha is required when use_thole = true.');
        end
        if ~isfield(sys, 'thole_a') || isempty(sys.thole_a)
            error('geom:build_active_row_cache_periodic:MissingTholeA', ...
                'sys.thole_a is required when use_thole = true.');
        end

        alphaFull = sys.site_alpha(:);
        if numel(alphaFull) ~= nSites
            error('geom:build_active_row_cache_periodic:BadSiteAlpha', ...
                'sys.site_alpha must have length equal to size(sys.site_pos,1).');
        end
        alphaAct = alphaFull(activeSites);

        a = sys.thole_a;
        if isscalar(a)
            tholeAAct = double(a);
        else
            a = a(:);
            if numel(a) ~= nSites
                error('geom:build_active_row_cache_periodic:BadTholeA', ...
                    'Non-scalar sys.thole_a must have length equal to size(sys.site_pos,1).');
            end
            tholeAAct = a(activeSites);
        end
    end

    t0 = tic;
    [row_ptr, col_idx, dr, r2_bare, r_bare, coeff_iso, coeff_dyad, nPairsUndirected] = ...
        mex_build_active_row_cache_periodic( ...
            spatial.frac, ...
            double(spatial.grid_shape), ...
            double(spatial.bin_head), ...
            double(spatial.bin_next), ...
            double(spatial.neighbor_offsets), ...
            H, ...
            double(rcut), ...
            double(alphaEwald), ...
            alphaAct, ...
            tholeAAct);
    profileTimes.query_pairs = toc(t0);

    profileTimes.undirected_geometry = 0.0;
    profileTimes.thole_undirected = 0.0;
    profileTimes.csr_counts = 0.0;
    profileTimes.csr_fill = 0.0;
    profileTimes.total = toc(tTotal);

    rowCache = struct();
    rowCache.mode = 'periodic_realspace_row';
    rowCache.activeSites = activeSites;
    rowCache.nActive = nActive;
    rowCache.row_ptr = row_ptr;
    rowCache.col_idx = col_idx;
    rowCache.dr = dr;
    rowCache.r2_bare = r2_bare;
    rowCache.r_bare = r_bare;
    rowCache.coeff_iso = coeff_iso;
    rowCache.coeff_dyad = coeff_dyad;
    rowCache.alpha = alphaEwald;
    rowCache.rcut = rcut;
    rowCache.use_thole = ~isempty(alphaAct);
    rowCache.nPairsUndirected = nPairsUndirected;
    rowCache.nEntriesDirected = numel(col_idx);
    rowCache.nInteractions = numel(col_idx);
    rowCache.isPeriodic = true;
    rowCache.lattice = H;
    rowCache.Lmin = Lmin;

    if doProfile
        rowCache.profileTimes = profileTimes;
        fprintf('build_active_row_cache_periodic timing summary:\n');
        fprintf('  build_spatial_index   : %.6f s\n', profileTimes.build_spatial_index);
        fprintf('  query_pairs           : %.6f s\n', profileTimes.query_pairs);
        fprintf('  undirected_geometry   : %.6f s\n', profileTimes.undirected_geometry);
        fprintf('  thole_undirected      : %.6f s\n', profileTimes.thole_undirected);
        fprintf('  csr_counts            : %.6f s\n', profileTimes.csr_counts);
        fprintf('  csr_fill              : %.6f s\n', profileTimes.csr_fill);
        fprintf('  total                 : %.6f s\n', profileTimes.total);
        fprintf('  nPairsUndirected      : %d\n', rowCache.nPairsUndirected);
        fprintf('  nEntriesDirected      : %d\n', rowCache.nEntriesDirected);
        fprintf('  Lmin                  : %.16e\n', Lmin);
        fprintf('  backend               : MEX\n');
    end
    return;
end

rowCache = local_build_active_row_cache_periodic_matlab(sys, problem, ewaldParams, opts, spatial, H, Lmin);
end

function rowCache = local_build_active_row_cache_periodic_matlab(sys, problem, ewaldParams, opts, spatial, H, Lmin)

doProfile = false;
if isfield(opts, 'profile') && ~isempty(opts.profile)
    doProfile = logical(opts.profile);
end

useThole = true;
if isfield(opts, 'use_thole') && ~isempty(opts.use_thole)
    useThole = logical(opts.use_thole);
end

tTotal = tic;

pos = sys.site_pos;
nSites = size(pos, 1);
activeSites = problem.activeSites(:);
nActive = numel(activeSites);

alphaEwald = ewaldParams.alpha;
rcut = ewaldParams.rcut;

profileTimes = struct();
profileTimes.build_spatial_index = 0.0;

queryOpts = struct();
queryOpts.return_r = true;
queryOpts.return_dr = true;
queryOpts.return_full_idx = false;

t0 = tic;
pairs = geom.query_pairs_within_cutoff(spatial, rcut, queryOpts);
profileTimes.query_pairs = toc(t0);

pairI = pairs.i(:);
pairJ = pairs.j(:);
drUndir = pairs.dr;
r2Undir = pairs.r2(:);
rUndir = pairs.r(:);

if any(r2Undir <= 0)
    error('geom:build_active_row_cache_periodic:NonPositiveDistance', ...
        'Encountered non-positive pair distance while building periodic row cache.');
end

% -------------------------------------------------------------------------
% DEDUPE UNDIRECTED PAIRS
% MATLAB periodic query returns the correct unique undirected pair set,
% but with duplicates on medium/large real systems. Canonicalize and keep
% first occurrence only.
% -------------------------------------------------------------------------
canonPairs = [min(pairI,pairJ), max(pairI,pairJ)];
[~, ia] = unique(canonPairs, 'rows', 'stable');

pairI = pairI(ia);
pairJ = pairJ(ia);
drUndir = drUndir(ia, :);
r2Undir = r2Undir(ia);
rUndir = rUndir(ia);

nPairs = numel(pairI);
nDir = 2 * nPairs;

t0 = tic;

alpha2 = alphaEwald^2;
twoAlphaOverSqrtPi = 2 * alphaEwald / sqrt(pi);

invR2 = 1 ./ r2Undir;
invR  = 1 ./ rUndir;
invR3 = invR .* invR2;
invR5 = invR3 .* invR2;
invR4 = invR2.^2;

erfcar = erfc(alphaEwald * rUndir);
expar2 = exp(-alpha2 * r2Undir);

B = erfcar .* invR3 + twoAlphaOverSqrtPi .* expar2 .* invR2;
C = 3 .* erfcar .* invR5 + ...
    twoAlphaOverSqrtPi .* (2 .* alpha2 .* invR2 + 3 .* invR4) .* expar2;

coeffIsoUndir = -B;
coeffDyadUndir = +C;

profileTimes.undirected_geometry = toc(t0);

t0 = tic;
haveThole = false;

if useThole
    if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
        error('geom:build_active_row_cache_periodic:MissingSiteAlpha', ...
            'sys.site_alpha is required when use_thole = true.');
    end
    if ~isfield(sys, 'thole_a') || isempty(sys.thole_a)
        error('geom:build_active_row_cache_periodic:MissingTholeA', ...
            'sys.thole_a is required when use_thole = true.');
    end

    alphaSite = sys.site_alpha(:);
    if numel(alphaSite) ~= nSites
        error('geom:build_active_row_cache_periodic:BadSiteAlpha', ...
            'sys.site_alpha must have length equal to size(sys.site_pos,1).');
    end

    alpha_i = alphaSite(activeSites(pairI));
    alpha_j = alphaSite(activeSites(pairJ));

    tf = local_thole_factors_vectorized(rUndir, alpha_i, alpha_j, sys.thole_a);

    coeffIsoUndir = coeffIsoUndir - tf.l3 .* invR3;
    coeffDyadUndir = coeffDyadUndir + 3 .* tf.l5 .* invR5;
    haveThole = true;
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
source_full_idx = zeros(nDir, 1);
drDir = zeros(nDir, 3);
r2Dir = zeros(nDir, 1);
rDir = zeros(nDir, 1);
coeffIsoDir = zeros(nDir, 1);
coeffDyadDir = zeros(nDir, 1);

nextPtr = row_ptr(1:end-1);

for p = 1:nPairs
    i = pairI(p);
    j = pairJ(p);

    k = nextPtr(i);
    col_idx(k) = j;
    source_full_idx(k) = activeSites(j);
    drDir(k, :) = drUndir(p, :);
    r2Dir(k) = r2Undir(p);
    rDir(k) = rUndir(p);
    coeffIsoDir(k) = coeffIsoUndir(p);
    coeffDyadDir(k) = coeffDyadUndir(p);
    nextPtr(i) = k + 1;

    k = nextPtr(j);
    col_idx(k) = i;
    source_full_idx(k) = activeSites(i);
    drDir(k, :) = -drUndir(p, :);
    r2Dir(k) = r2Undir(p);
    rDir(k) = rUndir(p);
    coeffIsoDir(k) = coeffIsoUndir(p);
    coeffDyadDir(k) = coeffDyadUndir(p);
    nextPtr(j) = k + 1;
end
profileTimes.csr_fill = toc(t0);

rowCache = struct();
rowCache.mode = 'periodic_realspace_row';
rowCache.activeSites = activeSites;
rowCache.nActive = nActive;
rowCache.row_ptr = row_ptr;
rowCache.col_idx = col_idx;
rowCache.source_full_idx = source_full_idx;
rowCache.dr = drDir;
rowCache.r2_bare = r2Dir;
rowCache.r_bare = rDir;
rowCache.coeff_iso = coeffIsoDir;
rowCache.coeff_dyad = coeffDyadDir;
rowCache.alpha = alphaEwald;
rowCache.rcut = rcut;
rowCache.use_thole = haveThole;
rowCache.nPairsUndirected = nPairs;
rowCache.nEntriesDirected = nDir;
rowCache.nInteractions = nDir;
rowCache.isPeriodic = true;
rowCache.lattice = H;
rowCache.Lmin = Lmin;

profileTimes.total = toc(tTotal);

if doProfile
    rowCache.profileTimes = profileTimes;
    fprintf('build_active_row_cache_periodic timing summary:\n');
    fprintf('  build_spatial_index   : %.6f s\n', profileTimes.build_spatial_index);
    fprintf('  query_pairs           : %.6f s\n', profileTimes.query_pairs);
    fprintf('  undirected_geometry   : %.6f s\n', profileTimes.undirected_geometry);
    fprintf('  thole_undirected      : %.6f s\n', profileTimes.thole_undirected);
    fprintf('  csr_counts            : %.6f s\n', profileTimes.csr_counts);
    fprintf('  csr_fill              : %.6f s\n', profileTimes.csr_fill);
    fprintf('  total                 : %.6f s\n', profileTimes.total);
    fprintf('  nPairsUndirected      : %d\n', nPairs);
    fprintf('  nEntriesDirected      : %d\n', nDir);
    fprintf('  Lmin                  : %.16e\n', Lmin);
    fprintf('  backend               : MATLAB\n');
end
end

function H = local_get_periodic_lattice(sys)
if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
    H = sys.super_lattice;
elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
    H = sys.lattice;
else
    error('geom:build_active_row_cache_periodic:MissingLattice', ...
        'Periodic row cache requires sys.super_lattice or sys.lattice.');
end
validateattributes(H, {'double'}, {'size',[3 3], 'finite', 'real'}, ...
    mfilename, 'periodic lattice');
end

function tf = local_thole_factors_vectorized(r, alpha_i, alpha_j, thole_a)
if any(alpha_i < 0) || any(alpha_j < 0)
    error('geom:build_active_row_cache_periodic:NegativeAlpha', ...
        'alpha_i and alpha_j must be nonnegative.');
end

if isscalar(thole_a)
    if thole_a < 0
        error('geom:build_active_row_cache_periodic:NegativeTholeA', ...
            'thole_a must be nonnegative.');
    end
    a = thole_a * ones(size(r));
else
    a = thole_a(:);
    if numel(a) ~= numel(r)
        error('geom:build_active_row_cache_periodic:VectorTholeASize', ...
            'Non-scalar thole_a must match the pair-array length in this helper.');
    end
end

r = r(:);
tf = struct();
tf.f3 = ones(size(r));
tf.f5 = ones(size(r));
tf.l3 = zeros(size(r));
tf.l5 = zeros(size(r));

mask = (a ~= 0) & (alpha_i ~= 0) & (alpha_j ~= 0);
if ~any(mask)
    return;
end

u = zeros(size(r));
u(mask) = r(mask) ./ (alpha_i(mask) .* alpha_j(mask)).^(1/6);
au3 = a(mask) .* u(mask).^3;

tf.f3(mask) = 1 - exp(-au3);
tf.f5(mask) = 1 - (1 + au3) .* exp(-au3);
tf.l3(mask) = tf.f3(mask) - 1;
tf.l5(mask) = tf.f5(mask) - 1;
end