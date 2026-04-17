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
%
% Notes
%   - Thole damping factors are evaluated using the bare intersite distance.
%   - If rcut is supplied, only pairs with bare separation |r_ij| <= rcut
%     are included.

    io.assert_atomic_units(sys);

    if nargin < 3 || isempty(dipoleParams)
        dipoleParams = struct();
    end

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('sys.site_pos is missing or empty.');
    end
    if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
        error('sys.site_alpha is missing or empty.');
    end
    if ~isfield(sys, 'thole_a') || isempty(sys.thole_a)
        error('sys.thole_a is missing or empty.');
    end

    pos = sys.site_pos;
    alpha = sys.site_alpha(:);
    a = sys.thole_a;

    nSites = size(pos, 1);

    if ~isequal(size(mu), [nSites, 3])
        error('mu must be N x 3.');
    end

    exclude_self = true;
    if isfield(dipoleParams, 'exclude_self') && ~isempty(dipoleParams.exclude_self)
        exclude_self = dipoleParams.exclude_self;
    end

    softening = 0.0;
    if isfield(dipoleParams, 'softening') && ~isempty(dipoleParams.softening)
        softening = dipoleParams.softening;
    end

    rcut2 = inf;
    if isfield(dipoleParams, 'rcut') && ~isempty(dipoleParams.rcut)
        rcut = dipoleParams.rcut;
        if ~isscalar(rcut) || rcut <= 0
            error('dipoleParams.rcut must be a positive scalar when provided.');
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
        error('target_mask must have length N.');
    end
    if numel(source_mask) ~= nSites
        error('source_mask must have length N.');
    end

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