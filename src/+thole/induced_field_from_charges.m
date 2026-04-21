function E = induced_field_from_charges(sys, fieldParams)
%INDUCED_FIELD_FROM_CHARGES Electric field at each site from point charges.
%
% Inputs
%   sys          canonical polarization-system struct in atomic units with fields:
%                .site_pos     N x 3 Cartesian positions
%                .site_charge  N x 1 charges
%
%   fieldParams  struct, optional fields:
%                .exclude_self           logical, default true
%                .softening              scalar, default 0
%                .rcut                   scalar cutoff in bohr, optional
%                .target_mask            N x 1 logical, optional
%                .source_mask            N x 1 logical, optional
%                .geom_cache             optional geometry cache
%                .use_thole_damping      logical, default false
%
% Output
%   E            N x 3 electric field at each target site
%
% Notes
%   Bare field formula:
%       E_i = sum_j q_j * r_ij / |r_ij|^3
%
%   If use_thole_damping = true, use damped charge->induced-dipole coupling:
%       E_i = sum_j q_j * f3_ij * r_ij / |r_ij|^3
%   where f3_ij is the Thole damping factor computed from
%   thole.thole_f3f5_factors(r, alpha_i, alpha_j, sys.thole_a).
%
%   With optional softening:
%       |r_ij|^2 -> |r_ij|^2 + softening^2
%
%   If rcut is supplied, only source-target pairs with bare separation
%       |r_ij| <= rcut
%   are included.
%
% Supported cache kinds
%   1) legacy unordered all-pairs cache from geom.build_nonperiodic_pair_cache
%   2) dedicated asymmetric target-source cache:
%        cache.kind = 'target_source_field_cache'

    if nargin < 2 || isempty(fieldParams)
        fieldParams = struct();
    end

    io.assert_atomic_units(sys);

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('calc:induced_field_from_charges:MissingSitePos', ...
            'sys.site_pos is missing or empty.');
    end

    if ~isfield(sys, 'site_charge') || isempty(sys.site_charge)
        error('calc:induced_field_from_charges:MissingSiteCharge', ...
            'sys.site_charge is missing or empty.');
    end

    pos = sys.site_pos;
    q = sys.site_charge(:);

    nSites = size(pos, 1);

    if numel(q) ~= nSites
        error('calc:induced_field_from_charges:BadChargeSize', ...
            'sys.site_charge must have length N.');
    end

    exclude_self = true;
    if isfield(fieldParams, 'exclude_self') && ~isempty(fieldParams.exclude_self)
        exclude_self = logical(fieldParams.exclude_self);
    end

    softening = 0.0;
    if isfield(fieldParams, 'softening') && ~isempty(fieldParams.softening)
        softening = fieldParams.softening;
    end
    if ~isscalar(softening) || ~isreal(softening) || ~isfinite(softening) || softening < 0
        error('calc:induced_field_from_charges:BadSoftening', ...
            'fieldParams.softening must be a nonnegative scalar.');
    end

    rcut = [];
    rcut2 = inf;
    if isfield(fieldParams, 'rcut') && ~isempty(fieldParams.rcut)
        rcut = fieldParams.rcut;
        if ~isscalar(rcut) || rcut <= 0 || ~isfinite(rcut)
            error('calc:induced_field_from_charges:BadRcut', ...
                'fieldParams.rcut must be a positive scalar when provided.');
        end
        rcut2 = rcut^2;
    end

    useThole = false;
    if isfield(fieldParams, 'use_thole_damping') && ~isempty(fieldParams.use_thole_damping)
        useThole = logical(fieldParams.use_thole_damping);
    end

    if useThole
        if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
            error('calc:induced_field_from_charges:MissingSiteAlpha', ...
                'sys.site_alpha is required when use_thole_damping = true.');
        end
        if ~isfield(sys, 'thole_a') || isempty(sys.thole_a)
            error('calc:induced_field_from_charges:MissingTholeA', ...
                'sys.thole_a is required when use_thole_damping = true.');
        end

        alphaSites = sys.site_alpha(:);
        if numel(alphaSites) ~= nSites
            error('calc:induced_field_from_charges:BadSiteAlpha', ...
                'sys.site_alpha must have length N.');
        end

        tholeA = sys.thole_a;
    else
        alphaSites = [];
        tholeA = [];
    end

    if isfield(fieldParams, 'target_mask') && ~isempty(fieldParams.target_mask)
        target_mask = logical(fieldParams.target_mask(:));
    elseif isfield(sys, 'target_mask') && ~isempty(sys.target_mask)
        target_mask = logical(sys.target_mask(:));
    else
        target_mask = true(nSites, 1);
    end

    if isfield(fieldParams, 'source_mask') && ~isempty(fieldParams.source_mask)
        source_mask = logical(fieldParams.source_mask(:));
    elseif isfield(sys, 'source_mask') && ~isempty(sys.source_mask)
        source_mask = logical(sys.source_mask(:));
    else
        source_mask = true(nSites, 1);
    end

    if numel(target_mask) ~= nSites
        error('calc:induced_field_from_charges:BadTargetMask', ...
            'target_mask must have length N.');
    end

    if numel(source_mask) ~= nSites
        error('calc:induced_field_from_charges:BadSourceMask', ...
            'source_mask must have length N.');
    end

    geomCache = [];
    if isfield(fieldParams, 'geom_cache') && ~isempty(fieldParams.geom_cache)
        geomCache = fieldParams.geom_cache;
    end

    if ~isempty(geomCache)
        if isfield(geomCache, 'kind') && strcmp(geomCache.kind, 'target_source_field_cache')
            E = local_apply_target_source_cache(q, alphaSites, tholeA, useThole, ...
                geomCache, target_mask, source_mask, exclude_self, softening, rcut);
        else
            E = local_apply_cached_legacy(q, alphaSites, tholeA, useThole, ...
                geomCache, target_mask, source_mask, exclude_self, softening, rcut);
        end
        return;
    end

    E = local_apply_direct_vectorized(pos, q, alphaSites, tholeA, useThole, ...
        target_mask, source_mask, exclude_self, softening, rcut2);
end

function E = local_apply_target_source_cache(q, alphaSites, tholeA, useThole, ...
    cache, target_mask, source_mask, exclude_self, softening, rcut)

    nSites = cache.nSites;
    E = zeros(nSites, 3);

    if ~isfield(cache, 'mode') || ~strcmp(cache.mode, 'nonperiodic')
        error('calc:induced_field_from_charges:BadTargetSourceCacheMode', ...
            'target_source_field_cache must have mode = ''nonperiodic''.');
    end

    if ~isempty(rcut)
        if ~isfield(cache, 'rcut') || isempty(cache.rcut) || isinf(cache.rcut) || cache.rcut + 1e-12 < rcut
            error('calc:induced_field_from_charges:CacheTooSmall', ...
                'fieldParams.geom_cache does not cover the requested rcut.');
        end
    end

    if isfield(cache, 'exclude_self') && logical(cache.exclude_self) ~= logical(exclude_self)
        error('calc:induced_field_from_charges:ExcludeSelfMismatch', ...
            'fieldParams.geom_cache.exclude_self does not match request.');
    end

    pair_i = cache.pair_target_full(:);
    pair_j = cache.pair_source_full(:);
    dr     = cache.dr;
    r2Bare = cache.r2_bare(:);

    keep = target_mask(pair_i) & source_mask(pair_j);

    if ~any(keep)
        return;
    end

    pair_i = pair_i(keep);
    pair_j = pair_j(keep);
    dr     = dr(keep, :);
    r2Bare = r2Bare(keep);

    qj = q(pair_j);
    keepq = (qj ~= 0);

    if ~any(keepq)
        return;
    end

    pair_i = pair_i(keepq);
    pair_j = pair_j(keepq);
    dr     = dr(keepq, :);
    r2Bare = r2Bare(keepq);
    qj     = qj(keepq);

    if softening == 0
        invR3 = cache.inv_r3_bare(keep);
        invR3 = invR3(keepq);
    else
        r2 = r2Bare + softening^2;
        invR = 1 ./ sqrt(r2);
        invR3 = invR ./ r2;
    end

    if useThole
        if softening == 0 && isfield(cache, 'thole_f3') && ~isempty(cache.thole_f3)
            f3 = cache.thole_f3(keep);
            f3 = f3(keepq);
        else
            nPairsLocal = numel(pair_i);
            f3 = zeros(nPairsLocal, 1);
            rBare = sqrt(r2Bare);
            for p = 1:nPairsLocal
                i = pair_i(p);
                j = pair_j(p);
                tf = thole.thole_f3f5_factors(rBare(p), alphaSites(i), alphaSites(j), tholeA);
                f3(p) = tf.f3;
            end
        end
        invR3 = invR3 .* f3;
    end

    coeff = qj .* invR3;
    contrib = (-dr) .* coeff;

    E(:, 1) = accumarray(pair_i, contrib(:, 1), [nSites, 1], @sum, 0);
    E(:, 2) = accumarray(pair_i, contrib(:, 2), [nSites, 1], @sum, 0);
    E(:, 3) = accumarray(pair_i, contrib(:, 3), [nSites, 1], @sum, 0);

    if ~exclude_self
        warning('calc:induced_field_from_charges:SelfNotCached', ...
            ['exclude_self = false requested with target-source cache, but self terms are ', ...
             'not represented unless explicitly included in the cache.']);
    end
end

function E = local_apply_cached_legacy(q, alphaSites, tholeA, useThole, ...
    cache, target_mask, source_mask, exclude_self, softening, rcut)

    nSites = cache.nSites;
    E = zeros(nSites, 3);

    if ~isempty(rcut)
        if ~isfield(cache, 'rcut') || isempty(cache.rcut) || isinf(cache.rcut) || cache.rcut + 1e-12 < rcut
            error('calc:induced_field_from_charges:CacheTooSmall', ...
                'fieldParams.geom_cache does not cover the requested rcut.');
        end
    end

    pair_i = cache.pair_i(:);
    pair_j = cache.pair_j(:);
    dr     = cache.dr;
    r2Bare = cache.r2_bare(:);

    if softening == 0
        invR3 = cache.inv_r3_bare(:);
    else
        r2 = r2Bare + softening^2;
        invR = 1 ./ sqrt(r2);
        invR3 = invR ./ r2;
    end

    if useThole
        if softening == 0 && isfield(cache, 'thole_f3') && ~isempty(cache.thole_f3)
            f3 = cache.thole_f3(:);
        else
            nPairsLocal = numel(pair_i);
            f3 = zeros(nPairsLocal, 1);
            rBare = sqrt(r2Bare);
            for p = 1:nPairsLocal
                i = pair_i(p);
                j = pair_j(p);
                tf = thole.thole_f3f5_factors(rBare(p), alphaSites(i), alphaSites(j), tholeA);
                f3(p) = tf.f3;
            end
        end
        invR3 = invR3 .* f3;
    end

    nPairs = numel(pair_i);

    for p = 1:nPairs
        i = pair_i(p);
        j = pair_j(p);
        rij = dr(p, :);  % r_j - r_i

        if target_mask(i) && source_mask(j)
            qj = q(j);
            if ~(exclude_self && i == j) && qj ~= 0
                E(i, :) = E(i, :) + qj * (-rij) * invR3(p);
            end
        end

        if target_mask(j) && source_mask(i)
            qi = q(i);
            if ~(exclude_self && i == j) && qi ~= 0
                E(j, :) = E(j, :) + qi * rij * invR3(p);
            end
        end
    end

    if ~exclude_self
        warning('calc:induced_field_from_charges:SelfNotCached', ...
            ['exclude_self = false requested with geom_cache, but cache only stores unordered ', ...
             'distinct pairs. Self terms are not included in cached path.']);
    end
end

function E = local_apply_direct_vectorized(pos, q, alphaSites, tholeA, useThole, ...
    target_mask, source_mask, exclude_self, softening, rcut2)

    nSites = size(pos, 1);
    E = zeros(nSites, 3);

    targetIdx = find(target_mask);
    if isempty(targetIdx)
        return;
    end

    % Only nonzero charges matter.
    sourceIdx = find(source_mask & (q ~= 0));
    if isempty(sourceIdx)
        return;
    end

    targetPos = pos(targetIdx, :);
    Ework = zeros(numel(targetIdx), 3);

    for b = 1:numel(sourceIdx)
        j = sourceIdx(b);
        qj = q(j);
        rj = pos(j, :);

        % r_ij = r_i - r_j for all targets
        rij = targetPos - rj;               % nTarget x 3
        r2_bare = sum(rij.^2, 2);          % nTarget x 1

        keep = true(numel(targetIdx), 1);

        if exclude_self
            keep = keep & (targetIdx ~= j);
        end

        if isfinite(rcut2)
            keep = keep & (r2_bare <= rcut2);
        end

        if ~any(keep)
            continue;
        end

        rij_keep = rij(keep, :);
        r2_bare_keep = r2_bare(keep);

        r2 = r2_bare_keep + softening^2;

        nonzero = (r2 > 0);
        if ~all(nonzero)
            rij_keep = rij_keep(nonzero, :);
            r2_bare_keep = r2_bare_keep(nonzero);
            keepIdx = find(keep);
            keepIdx = keepIdx(nonzero);
        else
            keepIdx = find(keep);
        end

        if isempty(keepIdx)
            continue;
        end

        invR = 1 ./ sqrt(r2_bare_keep + softening^2);
        invR3 = invR ./ (r2_bare_keep + softening^2);

        if useThole
            f3 = zeros(numel(keepIdx), 1);
            rBare = sqrt(r2_bare_keep);
            alpha_j = alphaSites(j);
            tgtFull = targetIdx(keepIdx);
            for k = 1:numel(keepIdx)
                i = tgtFull(k);
                tf = thole.thole_f3f5_factors(rBare(k), alphaSites(i), alpha_j, tholeA);
                f3(k) = tf.f3;
            end
            invR3 = invR3 .* f3;
        end

        Ework(keepIdx, :) = Ework(keepIdx, :) + qj .* rij_keep .* invR3;
    end

    E(targetIdx, :) = Ework;
end