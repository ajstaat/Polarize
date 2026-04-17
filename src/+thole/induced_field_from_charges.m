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
%                .geom_cache             optional bare-geometry cache from
%                                        geom.build_nonperiodic_pair_cache
%
% Output
%   E            N x 3 electric field at each target site
%
% Notes
%   Field formula used:
%       E_i = sum_j q_j * r_ij / |r_ij|^3
%   with optional softening:
%       |r_ij|^2 -> |r_ij|^2 + softening^2
%   If rcut is supplied, only source-target pairs with bare separation
%       |r_ij| <= rcut
%   are included.

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
        E = local_apply_cached(q, geomCache, target_mask, source_mask, exclude_self, softening, rcut);
        return;
    end

    E = local_apply_direct(pos, q, target_mask, source_mask, exclude_self, softening, rcut2);
end

function E = local_apply_cached(q, cache, target_mask, source_mask, exclude_self, softening, rcut)

    nSites = cache.nSites;
    E = zeros(nSites, 3);

    if ~isempty(rcut)
        if ~isfield(cache, 'rcut') || isempty(cache.rcut) || isinf(cache.rcut) || cache.rcut + 1e-12 < rcut
            error('calc:induced_field_from_charges:CacheTooSmall', ...
                'fieldParams.geom_cache does not cover the requested rcut.');
        end
    end

    if exclude_self && isequal(target_mask, source_mask)
        E = local_apply_cached_symmetric(E, q, cache, target_mask, softening);
    else
        E = local_apply_cached_masked(E, q, cache, target_mask, source_mask, softening, exclude_self);
    end
end

function E = local_apply_cached_symmetric(E, q, cache, mask, softening)
% Fast path for equal target/source masks with self excluded.
% Uses unordered pair symmetry: pair (i,j) updates both i and j.

    pair_i = cache.pair_i;
    pair_j = cache.pair_j;
    dr     = cache.dr;
    r2Bare = cache.r2_bare;

    usePair = mask(pair_i) & mask(pair_j);

    if ~any(usePair)
        return;
    end

    pair_i = pair_i(usePair);
    pair_j = pair_j(usePair);
    dr     = dr(usePair, :);
    r2Bare = r2Bare(usePair);

    if softening == 0
        invR3 = cache.inv_r3_bare(usePair);
    else
        r2 = r2Bare + softening^2;
        invR = 1 ./ sqrt(r2);
        invR3 = invR ./ r2;
    end

    nPairs = numel(pair_i);

    for p = 1:nPairs
        i = pair_i(p);
        j = pair_j(p);

        qi = q(i);
        qj = q(j);

        if qi == 0 && qj == 0
            continue;
        end

        rij = dr(p, :);  % r_j - r_i

        % Field at i from charge on j uses r_i - r_j = -rij
        if qj ~= 0
            E(i, :) = E(i, :) + qj * (-rij) * invR3(p);
        end

        % Field at j from charge on i uses r_j - r_i = +rij
        if qi ~= 0
            E(j, :) = E(j, :) + qi * rij * invR3(p);
        end
    end
end

function E = local_apply_cached_masked(E, q, cache, target_mask, source_mask, softening, exclude_self)
% General masked cached path.

    pair_i = cache.pair_i;
    pair_j = cache.pair_j;
    dr     = cache.dr;
    r2Bare = cache.r2_bare;

    if softening == 0
        invR3 = cache.inv_r3_bare;
    else
        r2 = r2Bare + softening^2;
        invR = 1 ./ sqrt(r2);
        invR3 = invR ./ r2;
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

function E = local_apply_direct(pos, q, target_mask, source_mask, exclude_self, softening, rcut2)

    nSites = size(pos, 1);
    E = zeros(nSites, 3);

    targetIdx = find(target_mask);
    sourceIdx = find(source_mask);

    for a = 1:numel(targetIdx)
        i = targetIdx(a);

        ri = pos(i, :);
        Ei = [0, 0, 0];

        for b = 1:numel(sourceIdx)
            j = sourceIdx(b);

            if exclude_self && i == j
                continue;
            end

            qj = q(j);
            if qj == 0
                continue;
            end

            rij = ri - pos(j, :);
            r2_bare = dot(rij, rij);

            if r2_bare > rcut2
                continue;
            end

            r2 = r2_bare + softening^2;

            if r2 == 0
                continue;
            end

            r3 = r2^(3/2);
            Ei = Ei + qj * rij / r3;
        end

        E(i, :) = Ei;
    end
end