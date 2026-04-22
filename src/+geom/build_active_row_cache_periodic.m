function rowCache = build_active_row_cache_periodic(sys, problem, ewaldParams, opts, spatial)
%BUILD_ACTIVE_ROW_CACHE_PERIODIC Build directed active-space row cache for
% periodic real-space dipole interactions using a periodic cell-list path.
%
% rowCache = geom.build_active_row_cache_periodic(sys, problem, ewaldParams)
% rowCache = geom.build_active_row_cache_periodic(sys, problem, ewaldParams, opts)
% rowCache = geom.build_active_row_cache_periodic(sys, problem, ewaldParams, opts, spatial)
%
% Fast path:
%   - uses periodic cell-list + MEX when available
%
% Fallback:
%   - uses existing direct periodic dipole-row builder
%
% Inputs
%   sys         canonical polarization system in atomic units
%   problem     struct with field .activeSites
%   ewaldParams struct with fields:
%                 .alpha
%                 .rcut
%
%   opts        optional struct with fields:
%                 .use_mex
%                 .profile
%                 .grid_shape    optional [gx gy gz]
%
%   spatial     optional prebuilt periodic spatial index on active-site positions
%
% Output
%   rowCache    struct with fields used by solve_scf_iterative_periodic_sor

narginchk(3, 5);
if nargin < 4 || isempty(opts)
    opts = struct();
end

io.assert_atomic_units(sys);

useMex = true;
if isfield(opts, 'use_mex') && ~isempty(opts.use_mex)
    useMex = logical(opts.use_mex);
end

doProfile = false;
if isfield(opts, 'profile') && ~isempty(opts.profile)
    doProfile = logical(opts.profile);
end

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
    error('geom:build_active_row_cache_periodic:MissingSitePos', ...
        'sys.site_pos is required.');
end
if ~isfield(problem, 'activeSites') || isempty(problem.activeSites)
    error('geom:build_active_row_cache_periodic:MissingActiveSites', ...
        'problem.activeSites is required.');
end
if ~isfield(ewaldParams, 'alpha') || isempty(ewaldParams.alpha)
    error('geom:build_active_row_cache_periodic:MissingAlpha', ...
        'ewaldParams.alpha is required.');
end
if ~isfield(ewaldParams, 'rcut') || isempty(ewaldParams.rcut)
    error('geom:build_active_row_cache_periodic:MissingRcut', ...
        'ewaldParams.rcut is required.');
end

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

alphaEwald = ewaldParams.alpha;
rcut = ewaldParams.rcut;

alphaFull = [];
if isfield(sys, 'site_alpha') && ~isempty(sys.site_alpha)
    alphaFull = sys.site_alpha(:);
    if numel(alphaFull) ~= nSites
        error('geom:build_active_row_cache_periodic:BadSiteAlpha', ...
            'sys.site_alpha must have length equal to size(sys.site_pos,1).');
    end
end
alphaAct = [];
if ~isempty(alphaFull)
    alphaAct = alphaFull(activeSites);
end

tholeAAct = [];
if isfield(sys, 'thole_a') && ~isempty(sys.thole_a)
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

H = local_get_direct_lattice(sys);

profileTimes = struct();
tTotal = tic;

% -------------------------------------------------------------------------
% Periodic spatial index on active sites
% -------------------------------------------------------------------------
if nargin < 5 || isempty(spatial)
    t0 = tic;
    spatial = local_build_periodic_spatial_index(posAct, H, rcut, opts);
    profileTimes.build_spatial_index = toc(t0);
else
    profileTimes.build_spatial_index = 0.0;
end

mexAvailable = (exist(['mex_build_active_row_cache_periodic.' mexext], 'file') == 3) || ...
               (exist('mex_build_active_row_cache_periodic', 'file') == 3);

useMexPath = useMex && mexAvailable && ...
    isstruct(spatial) && ...
    isfield(spatial, 'backend') && strcmp(spatial.backend, 'cell_list') && ...
    isfield(spatial, 'isPeriodic') && spatial.isPeriodic;

if useMexPath
    t0 = tic;
    [row_ptr, col_idx, source_full_idx_active, dr, r2_bare, r_bare, coeff_iso, coeff_dyad, ...
        nEntriesDirected, nCandidatesVisited, nCandidatesWithinCutoff] = ...
        mex_build_active_row_cache_periodic( ...
            spatial.pos, ...
            spatial.frac_pos_wrapped, ...
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
    profileTimes.total = toc(tTotal);

    rowCache = struct();
    rowCache.mode = 'periodic_realspace_dipole_row';
    rowCache.activeSites = activeSites;
    rowCache.targetSites = activeSites;
    rowCache.sourceSites = activeSites;
    rowCache.nActive = nActive;
    rowCache.nPolSites = nActive;
    rowCache.row_ptr = row_ptr;
    rowCache.col_idx = col_idx;
    rowCache.source_full_idx = activeSites(source_full_idx_active);
    rowCache.dr = dr;
    rowCache.r2_bare = r2_bare;
    rowCache.r_bare = r_bare;
    rowCache.coeff_iso = coeff_iso;
    rowCache.coeff_dyad = coeff_dyad;
    rowCache.rcut = rcut;
    rowCache.alpha = alphaEwald;
    rowCache.use_thole = ~isempty(tholeAAct);
    rowCache.nEntriesDirected = nEntriesDirected;
    rowCache.nInteractions = nEntriesDirected;
    rowCache.nCandidatesVisited = nCandidatesVisited;
    rowCache.nCandidatesWithinCutoff = nCandidatesWithinCutoff;

    if doProfile
        rowCache.profileTimes = profileTimes;
        fprintf('build_active_row_cache_periodic timing summary:\n');
        fprintf('  build_spatial_index   : %.6f s\n', profileTimes.build_spatial_index);
        fprintf('  query_pairs           : %.6f s\n', profileTimes.query_pairs);
        fprintf('  total                 : %.6f s\n', profileTimes.total);
        fprintf('  nEntriesDirected      : %d\n', rowCache.nEntriesDirected);
        fprintf('  nCandidatesVisited    : %d\n', rowCache.nCandidatesVisited);
        fprintf('  nWithinCutoff         : %d\n', rowCache.nCandidatesWithinCutoff);
        fprintf('  acceptance ratio      : %.6f\n', ...
            rowCache.nCandidatesWithinCutoff / max(rowCache.nCandidatesVisited, 1));
        fprintf('  avg candidates/row    : %.6f\n', ...
            rowCache.nCandidatesVisited / max(rowCache.nActive, 1));
        fprintf('  avg kept/row          : %.6f\n', ...
            rowCache.nEntriesDirected / max(rowCache.nActive, 1));
        fprintf('  backend               : MEX cell_list periodic\n');
    end
    return;
end

% -------------------------------------------------------------------------
% MATLAB fallback: use existing direct periodic dipole-row builder
% -------------------------------------------------------------------------
fallbackOpts = opts;
fallbackOpts.use_thole = ~isempty(tholeAAct);
fallbackOpts.profile = doProfile;

rowCache = geom.build_periodic_realspace_dipole_row_cache(sys, problem, ewaldParams, fallbackOpts);
rowCache.activeSites = activeSites;
rowCache.targetSites = activeSites;
rowCache.sourceSites = activeSites;
rowCache.nActive = nActive;
rowCache.nEntriesDirected = rowCache.nInteractions;
end

% =========================================================================
% Local periodic spatial index helper
%
% Metric-based stencil:
%   - use the lattice metric G = H' * H
%   - derive conservative image bounds from reciprocal vectors
%   - keep only raw bin offsets whose bin-box lower bound can be <= rcut
% =========================================================================
function spatial = local_build_periodic_spatial_index(posAct, H, rcut, opts)

nActive = size(posAct, 1);
frac = (H \ posAct.').';
frac = frac - floor(frac); % wrap into [0,1)

a = H(:,1); b = H(:,2); c = H(:,3);

% Choose a finer default grid than the failed coarse heuristic.
if isfield(opts, 'grid_shape') && ~isempty(opts.grid_shape)
    gridShape = round(opts.grid_shape(:).');
    if numel(gridShape) ~= 3 || any(gridShape < 1)
        error('geom:build_active_row_cache_periodic:BadGridShape', ...
            'opts.grid_shape must be a 1x3 vector of positive integers.');
    end
else
    % Aim for bin widths around rcut/2 along lattice-vector lengths.
    gridShape = max(1, ceil(2 * [norm(a) norm(b) norm(c)] / rcut));
end

gx = gridShape(1);
gy = gridShape(2);
gz = gridShape(3);
nBins = gx * gy * gz;

bx = min(floor(frac(:,1) * gx), gx - 1) + 1;
by = min(floor(frac(:,2) * gy), gy - 1) + 1;
bz = min(floor(frac(:,3) * gz), gz - 1) + 1;

binLin = bx + (by - 1) * gx + (bz - 1) * gx * gy;

bin_head = zeros(nBins, 1);
bin_next = zeros(nActive, 1);

for s = 1:nActive
    bLin = binLin(s);
    bin_next(s) = bin_head(bLin);
    bin_head(bLin) = s;
end

% Conservative image bounds from reciprocal basis.
% If b_i are the dual-lattice columns with b_i · a_j = delta_ij, then
% |n_i| <= rcut * ||b_i|| is a conservative bound whenever |H n| <= rcut.
Bdual = inv(H).';
nxmax = ceil(rcut * norm(Bdual(:,1))) + 1;
nymax = ceil(rcut * norm(Bdual(:,2))) + 1;
nzmax = ceil(rcut * norm(Bdual(:,3))) + 1;

neighbor_offsets = local_build_metric_neighbor_offsets(H, gridShape, ...
    nxmax, nymax, nzmax, rcut);

if isfield(opts, 'profile') && ~isempty(opts.profile) && logical(opts.profile)
    fprintf('build_active_row_cache_periodic spatial summary:\n');
    fprintf('  grid_shape           : [%d %d %d]\n', gridShape(1), gridShape(2), gridShape(3));
    fprintf('  image_bounds         : [%d %d %d]\n', nxmax, nymax, nzmax);
    fprintf('  nNeighborOffsets     : %d\n', size(neighbor_offsets, 1));
    fprintf('  nBins                : %d\n', nBins);
    fprintf('  nActive              : %d\n', nActive);
end

spatial = struct();
spatial.backend = 'cell_list';
spatial.isPeriodic = true;
spatial.pos = posAct;
spatial.frac_pos_wrapped = frac;
spatial.grid_shape = gridShape;
spatial.bin_head = bin_head;
spatial.bin_next = bin_next;
spatial.neighbor_offsets = neighbor_offsets;
end

function neighbor_offsets = local_build_metric_neighbor_offsets(H, gridShape, ...
    nxmax, nymax, nzmax, rcut)
% Build a conservative raw-bin stencil using a lattice-metric lower bound.
%
% Raw bin offsets encode both:
%   - home-cell bin displacement
%   - image shifts
%
% For raw offset d = [dx dy dz], the fractional center displacement is
%   c = [dx/gx, dy/gy, dz/gz]
%
% Points within the source/target bins can move relative to that center by
% e in [-hx, hx] x [-hy, hy] x [-hz, hz], where h = 1 ./ gridShape.
%
% We keep the offset if the lower bound
%   max(0, ||H c|| - rho)
% is <= rcut, where rho bounds max ||H e|| over e in the half-width box.

gx = gridShape(1);
gy = gridShape(2);
gz = gridShape(3);

G = H.' * H;
h = [1/gx, 1/gy, 1/gz];

% rho = max ||H e|| over e in [-h,h]^3
corners = [
    -1 -1 -1
    -1 -1  1
    -1  1 -1
    -1  1  1
     1 -1 -1
     1 -1  1
     1  1 -1
     1  1  1];
rho2 = 0.0;
for k = 1:size(corners,1)
    e = (corners(k,:) .* h).';
    rho2 = max(rho2, e.' * G * e);
end
rho = sqrt(rho2);

dxRange = -(nxmax * gx + (gx - 1)) : (nxmax * gx + (gx - 1));
dyRange = -(nymax * gy + (gy - 1)) : (nymax * gy + (gy - 1));
dzRange = -(nzmax * gz + (gz - 1)) : (nzmax * gz + (gz - 1));

offsets = zeros(numel(dxRange) * numel(dyRange) * numel(dzRange), 3);
nKeep = 0;

for dx = dxRange
    for dy = dyRange
        for dz = dzRange
            c = [dx / gx; dy / gy; dz / gz];
            dcenter2 = c.' * G * c;
            dcenter = sqrt(dcenter2);

            dmin = max(0.0, dcenter - rho);

            if dmin <= rcut
                nKeep = nKeep + 1;
                offsets(nKeep, :) = [dx dy dz];
            end
        end
    end
end

neighbor_offsets = offsets(1:nKeep, :);
end

function H = local_get_direct_lattice(sys)
if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
    H = sys.super_lattice;
elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
    H = sys.lattice;
else
    error('geom:build_active_row_cache_periodic:MissingLattice', ...
        'Missing direct lattice on system.');
end
end