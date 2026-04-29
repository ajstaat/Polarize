function rowCache = build_active_row_cache_periodic(sys, problem, ewaldParams, opts, spatial)
%BUILD_ACTIVE_ROW_CACHE_PERIODIC Build periodic real-space active row cache.
%
% rowCache = geom.build_active_row_cache_periodic(sys, problem, ewaldParams)
% rowCache = geom.build_active_row_cache_periodic(sys, problem, ewaldParams, opts)
% rowCache = geom.build_active_row_cache_periodic(sys, problem, ewaldParams, opts, spatial)
%
% Builds a directed active-space CSR row cache for the real-space part of
% periodic dipole Ewald interactions in the single-image regime.
%
% Project lattice convention
% --------------------------
% Polarize uses direct lattice vectors as ROWS:
%
%   r_cart = frac * H
%
% geom.get_lattice returns reciprocal vectors as COLUMNS:
%
%   H * G = 2*pi*I
%
% Public cache metadata follows this row-lattice convention. Any internal
% column/transpose operations are local implementation details.
%
% Inputs
%   sys           canonical polarization system with .site_pos and lattice
%   problem       struct with .activeSites
%   ewaldParams   struct with:
%                   .alpha
%                   .rcut
%   opts          optional struct:
%                   .use_mex     default true
%                   .profile     default false
%                   .use_thole   default true
%   spatial       optional prebuilt periodic spatial index for active sites
%
% Output
%   rowCache      struct with fields:
%                   .row_ptr
%                   .col_idx
%                   .dr
%                   .r2_bare
%                   .r_bare
%                   .coeff_iso
%                   .coeff_dyad
%
% Field action convention
% -----------------------
% For each directed row entry i <- j, the cached real-space contribution is
%
%   E_i += coeff_iso * mu_j + coeff_dyad * dot(mu_j, dr_ij) * dr_ij
%
% where dr_ij is the minimum-image vector from target i to source j.
%
% Notes
% -----
% This cache requires rcut < Lmin/2 so that a single minimum image is
% sufficient for real-space Ewald neighbor enumeration.

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

pos = sys.site_pos;
validateattributes(pos, {'numeric'}, {'2d', 'ncols', 3, 'real', 'finite'}, ...
    mfilename, 'sys.site_pos');
pos = double(pos);

nSites = size(pos, 1);

if ~isfield(problem, 'activeSites') || isempty(problem.activeSites)
    error('geom:build_active_row_cache_periodic:MissingActiveSites', ...
        'problem.activeSites is required and may not be empty.');
end

activeSites = problem.activeSites(:);
validateattributes(activeSites, {'numeric'}, ...
    {'vector', 'integer', 'positive', 'finite'}, ...
    mfilename, 'problem.activeSites');
activeSites = double(activeSites(:));

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

validateattributes(alphaEwald, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'ewaldParams.alpha');
validateattributes(rcut, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'ewaldParams.rcut');

alphaEwald = double(alphaEwald);
rcut = double(rcut);

lat = geom.get_lattice(sys);
H = lat.H;
G = lat.G;
V = lat.volume;

validateattributes(H, {'numeric'}, {'size', [3 3], 'real', 'finite'}, ...
    mfilename, 'periodic lattice');
validateattributes(G, {'numeric'}, {'size', [3 3], 'real', 'finite'}, ...
    mfilename, 'periodic reciprocal lattice');
validateattributes(V, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'periodic cell volume');

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

        alphaAct = double(alphaFull(activeSites));

        a = sys.thole_a;
        if isscalar(a)
            tholeAAct = double(a);
        else
            a = a(:);
            if numel(a) ~= nSites
                error('geom:build_active_row_cache_periodic:BadTholeA', ...
                    'Non-scalar sys.thole_a must have length equal to size(sys.site_pos,1).');
            end
            tholeAAct = double(a(activeSites));
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
    rowCache.source_full_idx = activeSites(col_idx);

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
    rowCache.lattice = lat;
    rowCache.H = H;
    rowCache.G = G;
    rowCache.V = V;
    rowCache.Lmin = Lmin;
    rowCache.lattice_convention = 'project_row_H_column_G_HG_2piI';
    if isfield(lat, 'convention')
        rowCache.convention = lat.convention;
    else
        rowCache.convention = 'project_row_H_column_G_HG_2piI';
    end

    rowCache.backend = 'MEX';

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

rowCache = local_build_active_row_cache_periodic_matlab( ...
    sys, problem, ewaldParams, opts, spatial, lat, Lmin);
end

function rowCache = local_build_active_row_cache_periodic_matlab( ...
    sys, problem, ewaldParams, opts, spatial, lat, Lmin)

doProfile = false;
if isfield(opts, 'profile') && ~isempty(opts.profile)
    doProfile = logical(opts.profile);
end

useThole = true;
if isfield(opts, 'use_thole') && ~isempty(opts.use_thole)
    useThole = logical(opts.use_thole);
end

tTotal = tic;

pos = double(sys.site_pos);
nSites = size(pos, 1);

activeSites = double(problem.activeSites(:));
nActive = numel(activeSites);

alphaEwald = double(ewaldParams.alpha);
rcut = double(ewaldParams.rcut);

H = lat.H;
G = lat.G;
V = lat.volume;

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
%
% The MATLAB periodic query can return duplicate undirected pairs on medium
% or large systems depending on the periodic cell-list stencil. Canonicalize
% local active-space pair labels and keep the first occurrence.
% -------------------------------------------------------------------------

canonPairs = [min(pairI, pairJ), max(pairI, pairJ)];
[~, ia] = unique(canonPairs, 'rows', 'stable');

pairI = pairI(ia);
pairJ = pairJ(ia);
drUndir = drUndir(ia, :);
r2Undir = r2Undir(ia);
rUndir = rUndir(ia);

nPairs = numel(pairI);
nDir = 2 * nPairs;

% -------------------------------------------------------------------------
% Periodic real-space dipole Ewald coefficients
%
% Field tensor form:
%
%   E_i += coeff_iso * mu_j + coeff_dyad * dot(mu_j, r_ij) * r_ij
%
% with
%
%   coeff_iso  = -B
%   coeff_dyad = +C
%
% and r_ij = source - target minimum-image vector.
% -------------------------------------------------------------------------

t0 = tic;

alpha2 = alphaEwald^2;
twoAlphaOverSqrtPi = 2 * alphaEwald / sqrt(pi);

invR2 = 1 ./ r2Undir;
invR = 1 ./ rUndir;
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

% -------------------------------------------------------------------------
% Optional short-range Thole correction.
%
% The real-space Ewald tensor above contains the undamped 1/r^3 singular
% short-range dipole tensor. Thole damping replaces the bare short-range
% tensor by the damped tensor, so we add the local correction:
%
%   T_damped - T_bare
%
% This is zero when use_thole = false.
% -------------------------------------------------------------------------

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

% -------------------------------------------------------------------------
% Build directed CSR row structure in local active-space indices.
% -------------------------------------------------------------------------

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
rowCache.lattice = lat;
rowCache.H = H;
rowCache.G = G;
rowCache.V = V;
rowCache.Lmin = Lmin;
rowCache.lattice_convention = 'project_row_H_column_G_HG_2piI';
if isfield(lat, 'convention')
    rowCache.convention = lat.convention;
else
    rowCache.convention = 'project_row_H_column_G_HG_2piI';
end

rowCache.backend = 'MATLAB';

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
alpha_i = alpha_i(:);
alpha_j = alpha_j(:);
a = a(:);

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

% Local correction factors relative to the bare tensor.
tf.l3(mask) = tf.f3(mask) - 1;
tf.l5(mask) = tf.f5(mask) - 1;
end