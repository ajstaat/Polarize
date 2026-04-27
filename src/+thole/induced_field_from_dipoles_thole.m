function Edip = induced_field_from_dipoles_thole(sys, mu, dipoleParams)
%INDUCED_FIELD_FROM_DIPOLES_THOLE Field at each site from induced dipoles
% with Thole damping consistent with thole_f3f5.
%
% Inputs
%   sys           canonical polarization-system struct in atomic units
%   mu            N x 3 induced dipoles
%   dipoleParams  optional struct with fields:
%                 .exclude_self   logical, default true
%                 .softening      scalar, default 0
%                 .rcut           scalar cutoff in bohr, optional
%                 .target_mask    N x 1 logical, optional
%                 .source_mask    N x 1 logical, optional
%                 .geom_cache     optional bare-geometry cache from
%                                 geom.build_nonperiodic_pair_cache
%
% Notes
%   - Thole damping factors are evaluated using the bare intersite distance.
%   - If rcut is supplied, only pairs with bare separation |r_ij| <= rcut
%     are included.
%   - When a geometry cache is supplied, the routine can avoid recomputing
%     pair geometry repeatedly.
%   - A symmetry-accelerated path is used when target/source masks are equal
%     and self-interactions are excluded.

    io.assert_atomic_units(sys);

    if nargin < 3 || isempty(dipoleParams)
        dipoleParams = struct();
    end

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('thole:induced_field_from_dipoles_thole:MissingSitePos', ...
            'sys.site_pos is missing or empty.');
    end
    if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
        error('thole:induced_field_from_dipoles_thole:MissingSiteAlpha', ...
            'sys.site_alpha is missing or empty.');
    end
    if ~isfield(sys, 'thole_a') || isempty(sys.thole_a)
        error('thole:induced_field_from_dipoles_thole:MissingTholeA', ...
            'sys.thole_a is missing or empty.');
    end

    pos = sys.site_pos;
    alpha = sys.site_alpha(:);
    a = sys.thole_a;

    nSites = size(pos, 1);

    if ~isequal(size(mu), [nSites, 3])
        error('thole:induced_field_from_dipoles_thole:BadMuSize', ...
            'mu must be N x 3.');
    end

    exclude_self = true;
    if isfield(dipoleParams, 'exclude_self') && ~isempty(dipoleParams.exclude_self)
        exclude_self = logical(dipoleParams.exclude_self);
    end

    softening = 0.0;
    if isfield(dipoleParams, 'softening') && ~isempty(dipoleParams.softening)
        softening = dipoleParams.softening;
    end
    if ~isscalar(softening) || ~isreal(softening) || ~isfinite(softening) || softening < 0
        error('thole:induced_field_from_dipoles_thole:BadSoftening', ...
            'dipoleParams.softening must be a nonnegative scalar.');
    end

    rcut = [];
    rcut2 = inf;
    if isfield(dipoleParams, 'rcut') && ~isempty(dipoleParams.rcut)
        rcut = dipoleParams.rcut;
        if ~isscalar(rcut) || rcut <= 0 || ~isfinite(rcut)
            error('thole:induced_field_from_dipoles_thole:BadRcut', ...
                'dipoleParams.rcut must be a positive scalar when provided.');
        end
        rcut2 = rcut^2;
    end

    if isfield(dipoleParams, 'target_mask') && ~isempty(dipoleParams.target_mask)
        target_mask = logical(dipoleParams.target_mask(:));
    else
        target_mask = true(nSites, 1);
    end

    if isfield(dipoleParams, 'source_mask') && ~isempty(dipoleParams.source_mask)
        source_mask = logical(dipoleParams.source_mask(:));
    else
        source_mask = true(nSites, 1);
    end

    if numel(target_mask) ~= nSites
        error('thole:induced_field_from_dipoles_thole:BadTargetMask', ...
            'target_mask must have length N.');
    end
    if numel(source_mask) ~= nSites
        error('thole:induced_field_from_dipoles_thole:BadSourceMask', ...
            'source_mask must have length N.');
    end

    geomCache = [];
    if isfield(dipoleParams, 'geom_cache') && ~isempty(dipoleParams.geom_cache)
        geomCache = dipoleParams.geom_cache;
    end

    if ~isempty(geomCache)
        Edip = local_apply_cached( ...
            mu, alpha, a, geomCache, ...
            target_mask, source_mask, ...
            exclude_self, softening, rcut);
        return;
    end

    Edip = local_apply_direct( ...
        pos, mu, alpha, a, ...
        target_mask, source_mask, ...
        exclude_self, softening, rcut2);
end

function Edip = local_apply_cached(mu, alpha, a, cache, target_mask, source_mask, ...
    exclude_self, softening, rcut)

    nSites = cache.nSites;
    Edip = zeros(nSites, 3);

    if ~isempty(rcut)
        if ~isfield(cache, 'rcut') || isempty(cache.rcut) || isinf(cache.rcut) || cache.rcut + 1e-12 < rcut
            error('thole:induced_field_from_dipoles_thole:CacheTooSmall', ...
                'dipoleParams.geom_cache does not cover the requested rcut.');
        end
    end

    if exclude_self && isequal(target_mask, source_mask)
        Edip = local_apply_cached_symmetric( ...
            Edip, mu, alpha, a, cache, target_mask, softening);
    else
        Edip = local_apply_cached_masked( ...
            Edip, mu, alpha, a, cache, target_mask, source_mask, softening, exclude_self);
    end
end

function Edip = local_apply_cached_symmetric(Edip, mu, alpha, a, cache, mask, softening)
% Fast path for equal target/source masks with self excluded.

    pair_i = cache.pair_i;
    pair_j = cache.pair_j;
    dr     = cache.dr;
    r2Bare = cache.r2_bare;
    rBare  = cache.r_bare;

    haveTholeFactors = isfield(cache, 'thole_f3') && isfield(cache, 'thole_f5');
    if haveTholeFactors
        tholeF3 = cache.thole_f3;
        tholeF5 = cache.thole_f5;
    end

    if softening == 0
        invR3 = cache.inv_r3_bare;
        invR5 = cache.inv_r5_bare;
    end

    usePair = mask(pair_i) & mask(pair_j);

    if ~any(usePair)
        return;
    end

    pair_i = pair_i(usePair);
    pair_j = pair_j(usePair);
    dr     = dr(usePair, :);
    r2Bare = r2Bare(usePair);
    rBare  = rBare(usePair);

    if haveTholeFactors
        tholeF3 = tholeF3(usePair);
        tholeF5 = tholeF5(usePair);
    end

    if softening == 0
        invR3 = invR3(usePair);
        invR5 = invR5(usePair);
    else
        r2 = r2Bare + softening^2;
        invR = 1 ./ sqrt(r2);
        invR3 = invR ./ r2;
        invR5 = invR3 ./ r2;
    end

    nPairs = numel(pair_i);

    for p = 1:nPairs
        i = pair_i(p);
        j = pair_j(p);

        muj = mu(j, :);
        mui = mu(i, :);

        if all(muj == 0) && all(mui == 0)
            continue;
        end

        rij = dr(p, :);  % r_j - r_i

        if haveTholeFactors
            f3 = tholeF3(p);
            f5 = tholeF5(p);
        else
            tf = thole.thole_f3f5_factors(rBare(p), alpha(i), alpha(j), a);
            f3 = tf.f3;
            f5 = tf.f5;
        end

        muDotR_ij = dot(muj, -rij);
        Edip(i, :) = Edip(i, :) ...
            + 3 * f5 * (-rij) * (muDotR_ij * invR5(p)) ...
            - f3 * muj * invR3(p);

        muDotR_ji = dot(mui, rij);
        Edip(j, :) = Edip(j, :) ...
            + 3 * f5 * rij * (muDotR_ji * invR5(p)) ...
            - f3 * mui * invR3(p);
    end
end

function Edip = local_apply_cached_masked(Edip, mu, alpha, a, cache, target_mask, ...
    source_mask, softening, exclude_self)
% General masked cached path.

    pair_i = cache.pair_i;
    pair_j = cache.pair_j;
    dr     = cache.dr;
    r2Bare = cache.r2_bare;
    rBare  = cache.r_bare;

    haveTholeFactors = isfield(cache, 'thole_f3') && isfield(cache, 'thole_f5');
    if haveTholeFactors
        tholeF3 = cache.thole_f3;
        tholeF5 = cache.thole_f5;
    end

    if softening == 0
        invR3 = cache.inv_r3_bare;
        invR5 = cache.inv_r5_bare;
    else
        r2 = r2Bare + softening^2;
        invR = 1 ./ sqrt(r2);
        invR3 = invR ./ r2;
        invR5 = invR3 ./ r2;
    end

    nPairs = numel(pair_i);

    for p = 1:nPairs
        i = pair_i(p);
        j = pair_j(p);
        rij = dr(p, :);

        if haveTholeFactors
            f3 = tholeF3(p);
            f5 = tholeF5(p);
        else
            tf = thole.thole_f3f5_factors(rBare(p), alpha(i), alpha(j), a);
            f3 = tf.f3;
            f5 = tf.f5;
        end

        if target_mask(i) && source_mask(j)
            muj = mu(j, :);
            if ~(exclude_self && i == j) && ~all(muj == 0)
                ri_minus_rj = -rij;
                muDotR = dot(muj, ri_minus_rj);
                Edip(i, :) = Edip(i, :) ...
                    + 3 * f5 * ri_minus_rj * (muDotR * invR5(p)) ...
                    - f3 * muj * invR3(p);
            end
        end

        if target_mask(j) && source_mask(i)
            mui = mu(i, :);
            if ~(exclude_self && i == j) && ~all(mui == 0)
                rj_minus_ri = rij;
                muDotR = dot(mui, rj_minus_ri);
                Edip(j, :) = Edip(j, :) ...
                    + 3 * f5 * rj_minus_ri * (muDotR * invR5(p)) ...
                    - f3 * mui * invR3(p);
            end
        end
    end

    if ~exclude_self
        warning('thole:induced_field_from_dipoles_thole:SelfNotCached', ...
            ['exclude_self = false requested with geom_cache, but cache only stores unordered ', ...
             'distinct pairs. Self terms are not included in cached path.']);
    end
end

function Edip = local_apply_direct(pos, mu, alpha, a, target_mask, source_mask, ...
    exclude_self, softening, rcut2)

    nSites = size(pos, 1);
    Edip = zeros(nSites, 3);

    targetIdx = find(target_mask);
    sourceIdx = find(source_mask);

    for aa = 1:numel(targetIdx)
        i = targetIdx(aa);
        ri = pos(i, :);
        Ei = [0, 0, 0];

        for bb = 1:numel(sourceIdx)
            j = sourceIdx(bb);

            if exclude_self && i == j
                continue;
            end

            muj = mu(j, :);
            if all(muj == 0)
                continue;
            end

            rij = ri - pos(j, :);
            r2_bare = dot(rij, rij);

            if r2_bare > rcut2
                continue;
            end

            if r2_bare == 0
                continue;
            end

            r2 = r2_bare + softening^2;
            r  = sqrt(r2);
            r3 = r2 * r;
            r5 = r2 * r3;

            muDotR = dot(muj, rij);

            r_phys = sqrt(r2_bare);
            tf = thole.thole_f3f5_factors(r_phys, alpha(i), alpha(j), a);
            Ei = Ei + 3 * tf.f5 * rij * (muDotR / r5) - tf.f3 * muj / r3;
        end

        Edip(i, :) = Ei;
    end
end