function cache = build_target_source_field_cache(sys, opts)
%BUILD_TARGET_SOURCE_FIELD_CACHE Build asymmetric target-source cache for
% nonperiodic external field from charges.
%
% cache = geom.build_target_source_field_cache(sys, opts)
%
% Inputs
%   sys  canonical polarization system in atomic units
%
%   opts struct with fields:
%     .target_mask   logical nSites x 1, required
%     .source_mask   logical nSites x 1, required
%     .exclude_self  logical scalar, default true
%     .rcut          scalar cutoff in bohr, default inf
%     .use_thole     logical scalar, default false
%     .softening     scalar, default 0.0
%
% Output
%   cache struct with fields:
%     .kind = 'target_source_field_cache'
%     .mode = 'nonperiodic'
%     .nSites
%     .nPairs
%     .pair_target_full
%     .pair_source_full
%     .dr = r_source - r_target
%     .r2_bare
%     .inv_r3_bare
%     .thole_f3 optional
%     .rcut
%     .exclude_self
%
% Notes
% - Only charged sources (sys.site_charge ~= 0) are retained as sources.
% - This is an asymmetric cache for external-field assembly, not a symmetric
%   dipole-dipole pair cache.
% - Neighbor generation is routed through the shared geometry tools by
%   querying a combined point set [targets; sources] and filtering to
%   cross-set pairs.

narginchk(2, 2);
io.assert_atomic_units(sys);

if ~isfield(opts, 'target_mask') || isempty(opts.target_mask)
    error('geom:build_target_source_field_cache:MissingTargetMask', ...
        'opts.target_mask is required.');
end
if ~isfield(opts, 'source_mask') || isempty(opts.source_mask)
    error('geom:build_target_source_field_cache:MissingSourceMask', ...
        'opts.source_mask is required.');
end
if ~isfield(sys, 'site_charge') || isempty(sys.site_charge)
    error('geom:build_target_source_field_cache:MissingSiteCharge', ...
        'sys.site_charge is required.');
end
if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
    error('geom:build_target_source_field_cache:MissingSitePos', ...
        'sys.site_pos is required.');
end

nSites = sys.n_sites;
pos = sys.site_pos;
q = sys.site_charge(:);

validateattributes(pos, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
    mfilename, 'sys.site_pos');
if size(pos, 1) ~= nSites
    error('geom:build_target_source_field_cache:SitePosSizeMismatch', ...
        'size(sys.site_pos,1) must equal sys.n_sites.');
end
if numel(q) ~= nSites
    error('geom:build_target_source_field_cache:SiteChargeSizeMismatch', ...
        'numel(sys.site_charge) must equal sys.n_sites.');
end

target_mask = logical(opts.target_mask(:));
source_mask = logical(opts.source_mask(:));
if numel(target_mask) ~= nSites || numel(source_mask) ~= nSites
    error('geom:build_target_source_field_cache:MaskSizeMismatch', ...
        'target_mask and source_mask must have length sys.n_sites.');
end

exclude_self = true;
if isfield(opts, 'exclude_self') && ~isempty(opts.exclude_self)
    exclude_self = logical(opts.exclude_self);
end

rcut = inf;
if isfield(opts, 'rcut') && ~isempty(opts.rcut)
    rcut = opts.rcut;
end
if ~(isscalar(rcut) && isreal(rcut) && isfinite(rcut) && rcut > 0) && ~isinf(rcut)
    error('geom:build_target_source_field_cache:BadRcut', ...
        'opts.rcut must be a positive scalar or Inf.');
end

useThole = false;
if isfield(opts, 'use_thole') && ~isempty(opts.use_thole)
    useThole = logical(opts.use_thole);
end

softening = 0.0;
if isfield(opts, 'softening') && ~isempty(opts.softening)
    softening = opts.softening;
end
if ~(isscalar(softening) && isreal(softening) && isfinite(softening) && softening >= 0)
    error('geom:build_target_source_field_cache:BadSoftening', ...
        'opts.softening must be a nonnegative finite scalar.');
end

targetIdx = find(target_mask);
% Only charged sources contribute to the external field.
sourceIdx = find(source_mask & (q ~= 0));

nTarget = numel(targetIdx);
nSource = numel(sourceIdx);

cache = struct();
cache.kind = 'target_source_field_cache';
cache.mode = 'nonperiodic';
cache.nSites = nSites;
cache.rcut = rcut;
cache.exclude_self = exclude_self;

if nTarget == 0 || nSource == 0
    cache.nPairs = 0;
    cache.pair_target_full = zeros(0, 1);
    cache.pair_source_full = zeros(0, 1);
    cache.dr = zeros(0, 3);
    cache.r2_bare = zeros(0, 1);
    cache.inv_r3_bare = zeros(0, 1);
    if useThole
        cache.thole_f3 = zeros(0, 1);
    end
    return;
end

targetPos = pos(targetIdx, :);
sourcePos = pos(sourceIdx, :);

% Combined point set: first targets, then sources.
comboPos = [targetPos; sourcePos];
nCombo = size(comboPos, 1);

if isfinite(rcut)
    spatial = geom.build_spatial_index(comboPos, struct( ...
        'isPeriodic', false, ...
        'method', 'auto', ...
        'cutoff', rcut));
else
    spatial = geom.build_spatial_index(comboPos, struct( ...
        'isPeriodic', false, ...
        'method', 'bruteforce'));
end

pairs = geom.query_pairs_within_cutoff(spatial, rcut, struct( ...
    'return_r', true, ...
    'return_dr', true, ...
    'return_full_idx', true));

pairI = pairs.i(:);
pairJ = pairs.j(:);
dr = pairs.dr;
r2_bare = pairs.r2(:);

isI_target = (pairI <= nTarget);
isJ_target = (pairJ <= nTarget);
isI_source = (pairI > nTarget);
isJ_source = (pairJ > nTarget);

% Keep only cross-set pairs: one target, one source.
keepCross = (isI_target & isJ_source) | (isI_source & isJ_target);

pairI = pairI(keepCross);
pairJ = pairJ(keepCross);
dr = dr(keepCross, :);
r2_bare = r2_bare(keepCross);

nCross = numel(pairI);

if nCross == 0
    cache.nPairs = 0;
    cache.pair_target_full = zeros(0, 1);
    cache.pair_source_full = zeros(0, 1);
    cache.dr = zeros(0, 3);
    cache.r2_bare = zeros(0, 1);
    cache.inv_r3_bare = zeros(0, 1);
    if useThole
        cache.thole_f3 = zeros(0, 1);
    end
    return;
end

pair_target_full = zeros(nCross, 1);
pair_source_full = zeros(nCross, 1);
dr_oriented = zeros(nCross, 3);

% Orient so that dr = r_source - r_target.
mask_target_source = (pairI <= nTarget) & (pairJ > nTarget);
if any(mask_target_source)
    iTS = pairI(mask_target_source);
    jTS = pairJ(mask_target_source);

    pair_target_full(mask_target_source) = targetIdx(iTS);
    pair_source_full(mask_target_source) = sourceIdx(jTS - nTarget);
    dr_oriented(mask_target_source, :) = dr(mask_target_source, :);
end

mask_source_target = (pairI > nTarget) & (pairJ <= nTarget);
if any(mask_source_target)
    iST = pairI(mask_source_target);
    jST = pairJ(mask_source_target);

    pair_target_full(mask_source_target) = targetIdx(jST);
    pair_source_full(mask_source_target) = sourceIdx(iST - nTarget);
    dr_oriented(mask_source_target, :) = -dr(mask_source_target, :);
end

if exclude_self
    keep = (pair_target_full ~= pair_source_full);
else
    keep = true(nCross, 1);
end

pair_target_full = pair_target_full(keep);
pair_source_full = pair_source_full(keep);
dr_oriented = dr_oriented(keep, :);
r2_bare = r2_bare(keep);

nPairs = numel(pair_target_full);

cache.nPairs = nPairs;

if nPairs == 0
    cache.pair_target_full = zeros(0, 1);
    cache.pair_source_full = zeros(0, 1);
    cache.dr = zeros(0, 3);
    cache.r2_bare = zeros(0, 1);
    cache.inv_r3_bare = zeros(0, 1);
    if useThole
        cache.thole_f3 = zeros(0, 1);
    end
    return;
end

cache.pair_target_full = pair_target_full;
cache.pair_source_full = pair_source_full;
cache.dr = dr_oriented;
cache.r2_bare = r2_bare;

r_bare = sqrt(r2_bare);
if softening == 0
    cache.inv_r3_bare = 1 ./ (r_bare.^3);
else
    % Keep the existing convention: softened denominators are handled
    % downstream; this cache stores bare geometry only.
    cache.inv_r3_bare = NaN(nPairs, 1);
end

if useThole
    if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
        error('geom:build_target_source_field_cache:MissingSiteAlpha', ...
            'sys.site_alpha is required when opts.use_thole is true.');
    end
    if ~isfield(sys, 'thole_a') || isempty(sys.thole_a)
        error('geom:build_target_source_field_cache:MissingTholeA', ...
            'sys.thole_a is required when opts.use_thole is true.');
    end

    alphaSites = sys.site_alpha(:);
    if numel(alphaSites) ~= nSites
        error('geom:build_target_source_field_cache:BadSiteAlpha', ...
            'numel(sys.site_alpha) must equal sys.n_sites.');
    end

    alpha_i = alphaSites(pair_target_full);
    alpha_j = alphaSites(pair_source_full);
    tf = thole.thole_f3f5_factors(r_bare, alpha_i, alpha_j, sys.thole_a);
    cache.thole_f3 = tf.f3;
end

end