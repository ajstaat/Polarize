function Edip = induced_field_from_dipoles_thole(sys, mu, dipoleParams)
%INDUCED_FIELD_FROM_DIPOLES_THOLE Field at each site from induced dipoles
% with Thole damping consistent with tholef3f5.

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

    target_mask = true(nSites, 1);
    if isfield(dipoleParams, 'target_mask') && ~isempty(dipoleParams.target_mask)
        target_mask = logical(dipoleParams.target_mask(:));
    end

    source_mask = true(nSites, 1);
    if isfield(dipoleParams, 'source_mask') && ~isempty(dipoleParams.source_mask)
        source_mask = logical(dipoleParams.source_mask(:));
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