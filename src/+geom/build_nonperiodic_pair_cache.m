function cache = build_nonperiodic_pair_cache(sys, opts, spatial)
%BUILD_NONPERIODIC_PAIR_CACHE Build bare geometry cache for nonperiodic pairs.
%
% cache = geom.build_nonperiodic_pair_cache(sys)
% cache = geom.build_nonperiodic_pair_cache(sys, opts)
% cache = geom.build_nonperiodic_pair_cache(sys, opts, spatial)
%
% Inputs
%   sys : struct with field
%         .site_pos   N x 3 Cartesian site coordinates
%
%   opts : optional struct with fields
%          .rcut      positive scalar cutoff in bohr, default = Inf
%          .site_mask N x 1 logical or numeric mask selecting sites to include
%          .use_mex   logical, default true
%          .profile   logical, default false
%
%   spatial : optional struct from geom.build_spatial_index(sys.site_pos, ...)
%
% Output
%   cache : struct with full-space unordered pair geometry
%           .nSites
%           .site_mask
%           .site_idx
%           .pair_i
%           .pair_j
%           .dr
%           .r2_bare
%           .r_bare
%           .inv_r_bare
%           .inv_r3_bare
%           .inv_r5_bare
%           .thole_f3 (if sys.site_alpha and sys.thole_a are present)
%           .thole_f5 (if sys.site_alpha and sys.thole_a are present)
%           .rcut
%           .nPairs
%           .profileTimes (if opts.profile = true)
%
% Notes
% - pair_i and pair_j are FULL site indices into sys.site_pos.
% - This cache stores bare geometry and optional Thole pair factors.
% - Softened denominators remain runtime-dependent elsewhere.
% - Uses MEX automatically when available for finite-rcut nonperiodic
%   cell-list builds, otherwise falls back to the preserved MATLAB version.

narginchk(1, 3);
if nargin < 2 || isempty(opts)
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

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
    error('geom:build_nonperiodic_pair_cache:MissingSitePos', ...
        'sys.site_pos is required and may not be empty.');
end

pos = sys.site_pos;
validateattributes(pos, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
    mfilename, 'sys.site_pos');

nSites = size(pos, 1);

rcut = local_get_opt(opts, 'rcut', inf);
if ~(isscalar(rcut) && isreal(rcut) && isfinite(rcut) && rcut > 0) && ~isinf(rcut)
    error('geom:build_nonperiodic_pair_cache:BadRcut', ...
        'opts.rcut must be a positive scalar or Inf.');
end

siteMask = local_get_opt(opts, 'site_mask', true(nSites, 1));
siteMask = logical(siteMask(:));
if numel(siteMask) ~= nSites
    error('geom:build_nonperiodic_pair_cache:BadSiteMask', ...
        'opts.site_mask must have length equal to size(sys.site_pos,1).');
end

profileTimes = struct();
tTotal = tic;

% -------------------------------------------------------------------------
% Spatial index
% -------------------------------------------------------------------------
if nargin < 3 || isempty(spatial)
    spatialOpts = struct();
    spatialOpts.isPeriodic = false;

    if isfinite(rcut)
        spatialOpts.method = 'auto';
        spatialOpts.cutoff = rcut;
    else
        spatialOpts.method = 'bruteforce';
    end

    t0 = tic;
    spatial = geom.build_spatial_index(pos, spatialOpts);
    profileTimes.build_spatial_index = toc(t0);
else
    profileTimes.build_spatial_index = 0.0;
end

siteIdx = find(siteMask);

mexAvailable = (exist(['mex_build_nonperiodic_pair_cache.' mexext], 'file') == 3) || ...
               (exist('mex_build_nonperiodic_pair_cache', 'file') == 3);

useMexPath = useMex && mexAvailable && ...
    isfinite(rcut) && ...
    isstruct(spatial) && ...
    isfield(spatial, 'backend') && strcmp(spatial.backend, 'cell_list') && ...
    isfield(spatial, 'isPeriodic') && ~spatial.isPeriodic;

if useMexPath
    alpha = [];
    tholeA = [];

    if isfield(sys, 'site_alpha') && ~isempty(sys.site_alpha) && ...
       isfield(sys, 'thole_a') && ~isempty(sys.thole_a)

        alpha = sys.site_alpha(:);
        if numel(alpha) ~= nSites
            error('geom:build_nonperiodic_pair_cache:BadSiteAlpha', ...
                'sys.site_alpha must have length equal to size(sys.site_pos,1).');
        end
        tholeA = sys.thole_a;
    end

    t0 = tic;
    [pair_i, pair_j, dr, r2Bare, rBare, invRBare, invR3Bare, invR5Bare, thole_f3, thole_f5, nPairs] = ...
        mex_build_nonperiodic_pair_cache( ...
            spatial.pos, ...
            double(siteMask), ...
            double(spatial.grid_shape), ...
            double(spatial.bin_head), ...
            double(spatial.bin_next), ...
            double(spatial.neighbor_offsets), ...
            double(rcut), ...
            alpha, ...
            tholeA);
    profileTimes.query_pairs = toc(t0);

    profileTimes.geometry = 0.0;
    profileTimes.thole = 0.0;
    profileTimes.total = toc(tTotal);

    cache = struct();
    cache.nSites = nSites;
    cache.site_mask = siteMask;
    cache.site_idx = siteIdx;
    cache.pair_i = pair_i;
    cache.pair_j = pair_j;
    cache.dr = dr;
    cache.r2_bare = r2Bare;
    cache.r_bare = rBare;
    cache.inv_r_bare = invRBare;
    cache.inv_r3_bare = invR3Bare;
    cache.inv_r5_bare = invR5Bare;
    cache.rcut = rcut;
    cache.nPairs = nPairs;

    if ~isempty(thole_f3)
        cache.thole_f3 = thole_f3;
        cache.thole_f5 = thole_f5;
    end

    if doProfile
        cache.profileTimes = profileTimes;
        fprintf('build_nonperiodic_pair_cache timing summary:\n');
        fprintf('  build_spatial_index   : %.6f s\n', profileTimes.build_spatial_index);
        fprintf('  query_pairs           : %.6f s\n', profileTimes.query_pairs);
        fprintf('  geometry              : %.6f s\n', profileTimes.geometry);
        fprintf('  thole                 : %.6f s\n', profileTimes.thole);
        fprintf('  total                 : %.6f s\n', profileTimes.total);
        fprintf('  nPairs                : %d\n', cache.nPairs);
        fprintf('  backend               : MEX\n');
    end
    return;
end

% -------------------------------------------------------------------------
% MATLAB fallback: preserve pasted implementation structure
% -------------------------------------------------------------------------
cache = local_build_nonperiodic_pair_cache_matlab(sys, opts, spatial, ...
    nSites, rcut, siteMask, siteIdx, doProfile, tTotal, profileTimes);

end

% =========================================================================
% MATLAB fallback (preserved from pasted version, with light profiling wrap)
% =========================================================================
function cache = local_build_nonperiodic_pair_cache_matlab(sys, opts, spatial, ...
    nSites, rcut, siteMask, siteIdx, doProfile, tTotal, profileTimes)

% -------------------------------------------------------------------------
% Direct finite-rcut cell-list path
% -------------------------------------------------------------------------
useDirectCellList = isfinite(rcut) && ...
    isstruct(spatial) && ...
    isfield(spatial, 'backend') && strcmp(spatial.backend, 'cell_list') && ...
    isfield(spatial, 'isPeriodic') && ~spatial.isPeriodic;

if useDirectCellList
    t0 = tic;
    [pair_i, pair_j, dr, r2Bare, rBare] = local_build_pairs_from_cell_list(spatial, siteMask, rcut);
    profileTimes.query_pairs = toc(t0);
else
    % ---------------------------------------------------------------------
    % Fallback path: brute-force query semantics
    % ---------------------------------------------------------------------
    queryOpts = struct();
    queryOpts.subset_idx = siteIdx;
    queryOpts.return_r = true;
    queryOpts.return_dr = true;
    queryOpts.return_full_idx = true;

    t0 = tic;
    pairs = geom.query_pairs_within_cutoff(spatial, rcut, queryOpts);
    profileTimes.query_pairs = toc(t0);

    pair_i = pairs.i(:);
    pair_j = pairs.j(:);
    dr = pairs.dr;
    r2Bare = pairs.r2(:);
    rBare = pairs.r(:);
end

t0 = tic;
if any(r2Bare <= 0)
    error('geom:build_nonperiodic_pair_cache:NonPositiveDistance', ...
        'Encountered non-positive pair distance while building cache.');
end

invRBare = 1 ./ rBare;
invR3Bare = invRBare ./ r2Bare;
invR5Bare = invR3Bare ./ r2Bare;
profileTimes.geometry = toc(t0);

cache = struct();
cache.nSites = nSites;
cache.site_mask = siteMask;
cache.site_idx = siteIdx;
cache.pair_i = pair_i;
cache.pair_j = pair_j;
cache.dr = dr;
cache.r2_bare = r2Bare;
cache.r_bare = rBare;
cache.inv_r_bare = invRBare;
cache.inv_r3_bare = invR3Bare;
cache.inv_r5_bare = invR5Bare;

t0 = tic;
if isfield(sys, 'site_alpha') && ~isempty(sys.site_alpha) && ...
   isfield(sys, 'thole_a') && ~isempty(sys.thole_a)

    alpha = sys.site_alpha(:);
    if numel(alpha) ~= nSites
        error('geom:build_nonperiodic_pair_cache:BadSiteAlpha', ...
            'sys.site_alpha must have length equal to size(sys.site_pos,1).');
    end

    tf = thole.thole_f3f5_factors(rBare, alpha(pair_i), alpha(pair_j), sys.thole_a);
    cache.thole_f3 = tf.f3;
    cache.thole_f5 = tf.f5;
end
profileTimes.thole = toc(t0);

cache.rcut = rcut;
cache.nPairs = numel(pair_i);

profileTimes.total = toc(tTotal);

if doProfile
    cache.profileTimes = profileTimes;
    fprintf('build_nonperiodic_pair_cache timing summary:\n');
    fprintf('  build_spatial_index   : %.6f s\n', profileTimes.build_spatial_index);
    fprintf('  query_pairs           : %.6f s\n', profileTimes.query_pairs);
    fprintf('  geometry              : %.6f s\n', profileTimes.geometry);
    fprintf('  thole                 : %.6f s\n', profileTimes.thole);
    fprintf('  total                 : %.6f s\n', profileTimes.total);
    fprintf('  nPairs                : %d\n', cache.nPairs);
    fprintf('  backend               : MATLAB\n');
end

end

% =========================================================================
% Preserved direct cell-list traversal from pasted fallback
% =========================================================================
function [pair_i, pair_j, dr_all, r2_all, r_all] = local_build_pairs_from_cell_list(spatial, siteMask, rcut)
cutoff2 = rcut^2;
gridShape = spatial.grid_shape;
offsets = spatial.neighbor_offsets;
nOffsets = size(offsets, 1);

% Pass 1: count pairs exactly
nPairs = 0;

for binLin = 1:spatial.n_bins
    iHead = spatial.bin_head(binLin);
    if iHead == 0
        continue;
    end

    [ix, iy, iz] = local_unlinearize(binLin, gridShape);

    for oo = 1:nOffsets
        off = offsets(oo, :);
        jx = ix + off(1);
        jy = iy + off(2);
        jz = iz + off(3);

        if jx < 1 || jx > gridShape(1) || ...
           jy < 1 || jy > gridShape(2) || ...
           jz < 1 || jz > gridShape(3)
            continue;
        end

        jBin = local_linearize([jx jy jz], gridShape);
        jHead = spatial.bin_head(jBin);
        if jHead == 0
            continue;
        end

        if jBin == binLin
            i = iHead;
            while i ~= 0
                if siteMask(i)
                    j = spatial.bin_next(i);
                    while j ~= 0
                        if siteMask(j)
                            dr = spatial.pos(j, :) - spatial.pos(i, :);
                            r2 = dot(dr, dr);
                            if r2 <= cutoff2
                                nPairs = nPairs + 1;
                            end
                        end
                        j = spatial.bin_next(j);
                    end
                end
                i = spatial.bin_next(i);
            end
        else
            i = iHead;
            while i ~= 0
                if siteMask(i)
                    j = jHead;
                    while j ~= 0
                        if siteMask(j)
                            dr = spatial.pos(j, :) - spatial.pos(i, :);
                            r2 = dot(dr, dr);
                            if r2 <= cutoff2
                                nPairs = nPairs + 1;
                            end
                        end
                        j = spatial.bin_next(j);
                    end
                end
                i = spatial.bin_next(i);
            end
        end
    end
end

pair_i = zeros(nPairs, 1);
pair_j = zeros(nPairs, 1);
dr_all = zeros(nPairs, 3);
r2_all = zeros(nPairs, 1);
r_all = zeros(nPairs, 1);

% Pass 2: fill
k = 0;

for binLin = 1:spatial.n_bins
    iHead = spatial.bin_head(binLin);
    if iHead == 0
        continue;
    end

    [ix, iy, iz] = local_unlinearize(binLin, gridShape);

    for oo = 1:nOffsets
        off = offsets(oo, :);
        jx = ix + off(1);
        jy = iy + off(2);
        jz = iz + off(3);

        if jx < 1 || jx > gridShape(1) || ...
           jy < 1 || jy > gridShape(2) || ...
           jz < 1 || jz > gridShape(3)
            continue;
        end

        jBin = local_linearize([jx jy jz], gridShape);
        jHead = spatial.bin_head(jBin);
        if jHead == 0
            continue;
        end

        if jBin == binLin
            i = iHead;
            while i ~= 0
                if siteMask(i)
                    j = spatial.bin_next(i);
                    while j ~= 0
                        if siteMask(j)
                            dr = spatial.pos(j, :) - spatial.pos(i, :);
                            r2 = dot(dr, dr);
                            if r2 <= cutoff2
                                k = k + 1;
                                pair_i(k) = i;
                                pair_j(k) = j;
                                dr_all(k, :) = dr;
                                r2_all(k) = r2;
                                r_all(k) = sqrt(r2);
                            end
                        end
                        j = spatial.bin_next(j);
                    end
                end
                i = spatial.bin_next(i);
            end
        else
            i = iHead;
            while i ~= 0
                if siteMask(i)
                    j = jHead;
                    while j ~= 0
                        if siteMask(j)
                            dr = spatial.pos(j, :) - spatial.pos(i, :);
                            r2 = dot(dr, dr);
                            if r2 <= cutoff2
                                k = k + 1;
                                pair_i(k) = i;
                                pair_j(k) = j;
                                dr_all(k, :) = dr;
                                r2_all(k) = r2;
                                r_all(k) = sqrt(r2);
                            end
                        end
                        j = spatial.bin_next(j);
                    end
                end
                i = spatial.bin_next(i);
            end
        end
    end
end

if k ~= nPairs
    error('geom:build_nonperiodic_pair_cache:PairCountMismatch', ...
        'Internal pair-count mismatch in direct cell-list traversal.');
end
end

% =========================================================================
% Helpers
% =========================================================================
function lin = local_linearize(ijk, gridShape)
lin = ijk(1) + (ijk(2)-1)*gridShape(1) + (ijk(3)-1)*gridShape(1)*gridShape(2);
end

function [ix, iy, iz] = local_unlinearize(lin, gridShape)
lin0 = lin - 1;
nxy = gridShape(1) * gridShape(2);
iz = floor(lin0 / nxy) + 1;
rem1 = lin0 - (iz - 1) * nxy;
iy = floor(rem1 / gridShape(1)) + 1;
ix = rem1 - (iy - 1) * gridShape(1) + 1;
end

function value = local_get_opt(s, name, defaultValue)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end