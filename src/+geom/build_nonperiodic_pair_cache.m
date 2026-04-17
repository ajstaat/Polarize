function cache = build_nonperiodic_pair_cache(sys, opts, spatial)
%BUILD_NONPERIODIC_PAIR_CACHE Build bare geometry cache for nonperiodic pairs.
%
% cache = geom.build_nonperiodic_pair_cache(sys)
% cache = geom.build_nonperiodic_pair_cache(sys, opts)
% cache = geom.build_nonperiodic_pair_cache(sys, opts, spatial)
%
% Inputs
%   sys : struct with field
%       .site_pos   N x 3 Cartesian site coordinates
%
%   opts : optional struct with fields
%       .rcut        positive scalar cutoff in bohr, default = Inf
%       .site_mask   N x 1 logical or numeric mask selecting sites to include
%
%   spatial : optional struct from geom.build_spatial_index(sys.site_pos, ...)
%
% Output
%   cache : struct with full-space unordered pair geometry
%       .nSites
%       .site_mask
%       .site_idx
%       .pair_i
%       .pair_j
%       .dr
%       .r2_bare
%       .r_bare
%       .inv_r_bare
%       .inv_r3_bare
%       .inv_r5_bare
%       .thole_f3        (if sys.site_alpha and sys.thole_a are present)
%       .thole_f5        (if sys.site_alpha and sys.thole_a are present)
%       .rcut
%       .nPairs
%
% Notes
%   - pair_i and pair_j are FULL site indices into sys.site_pos.
%   - This cache stores bare geometry and optional Thole pair factors.
%   - Softened denominators should remain runtime-dependent in +thole.

    narginchk(1, 3);

    if nargin < 2 || isempty(opts)
        opts = struct();
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

    if nargin < 3 || isempty(spatial)
        spatial = geom.build_spatial_index(pos, struct('isPeriodic', false));
    end

    siteIdx = find(siteMask);

    queryOpts = struct();
    queryOpts.subset_idx = siteIdx;
    queryOpts.return_r = true;
    queryOpts.return_dr = true;
    queryOpts.return_full_idx = true;

    pairs = geom.query_pairs_within_cutoff(spatial, rcut, queryOpts);

    r2Bare = pairs.r2;
    rBare  = pairs.r;

    if any(r2Bare <= 0)
        error('geom:build_nonperiodic_pair_cache:NonPositiveDistance', ...
            'Encountered non-positive pair distance while building cache.');
    end

    invRBare  = 1 ./ rBare;
    invR3Bare = invRBare ./ r2Bare;
    invR5Bare = invR3Bare ./ r2Bare;

    cache = struct();
    cache.nSites = nSites;
    cache.site_mask = siteMask;
    cache.site_idx = siteIdx;

    cache.pair_i = pairs.i;
    cache.pair_j = pairs.j;
    cache.dr = pairs.dr;

    cache.r2_bare = r2Bare;
    cache.r_bare = rBare;
    cache.inv_r_bare = invRBare;
    cache.inv_r3_bare = invR3Bare;
    cache.inv_r5_bare = invR5Bare;

    % Optional precomputed Thole factors for repeated field applications.
    if isfield(sys, 'site_alpha') && ~isempty(sys.site_alpha) && ...
       isfield(sys, 'thole_a') && ~isempty(sys.thole_a)

        alpha = sys.site_alpha(:);
        if numel(alpha) ~= nSites
            error('geom:build_nonperiodic_pair_cache:BadSiteAlpha', ...
                'sys.site_alpha must have length equal to size(sys.site_pos,1).');
        end

        a = sys.thole_a;
        nPairs = numel(pairs.i);
        tholeF3 = zeros(nPairs, 1);
        tholeF5 = zeros(nPairs, 1);

        for p = 1:nPairs
            i = pairs.i(p);
            j = pairs.j(p);
            tf = thole.thole_f3f5_factors(rBare(p), alpha(i), alpha(j), a);
            tholeF3(p) = tf.f3;
            tholeF5(p) = tf.f5;
        end

        cache.thole_f3 = tholeF3;
        cache.thole_f5 = tholeF5;
    end

    cache.rcut = rcut;
    cache.nPairs = numel(pairs.i);
end

function value = local_get_opt(s, name, defaultValue)
    if isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = defaultValue;
    end
end