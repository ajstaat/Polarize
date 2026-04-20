function coeffCache = build_periodic_realspace_dipole_coeff_cache(sys, rowGeomCache, ewaldParams, scfParams)
%BUILD_PERIODIC_REALSPACE_DIPOLE_COEFF_CACHE Build dipole-specific real-space coefficients
% on top of shared row geometry cache.
%
% coeffCache = geom.build_periodic_realspace_dipole_coeff_cache(sys, rowGeomCache, ewaldParams, scfParams)
%
% Output fields:
%   .coeff_iso
%   .coeff_dyad
%   .alpha
%   .rcut
%   .use_thole

    narginchk(4, 4);

    if nargin < 4 || isempty(scfParams)
        scfParams = struct();
    end

    use_thole = true;
    if isfield(scfParams, 'use_thole') && ~isempty(scfParams.use_thole)
        use_thole = logical(scfParams.use_thole);
    end

    alpha = ewaldParams.alpha;
    validateattributes(alpha, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'ewaldParams.alpha');

    nSites = rowGeomCache.nSites;
    if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
        error('geom:build_periodic_realspace_dipole_coeff_cache:MissingSiteAlpha', ...
            'sys.site_alpha is required.');
    end
    alpha_site = sys.site_alpha(:);
    if numel(alpha_site) ~= nSites
        error('geom:build_periodic_realspace_dipole_coeff_cache:BadSiteAlpha', ...
            'sys.site_alpha must have length equal to nSites.');
    end

    if use_thole
        if ~isfield(sys, 'thole_a') || isempty(sys.thole_a)
            error('geom:build_periodic_realspace_dipole_coeff_cache:MissingTholeA', ...
                'sys.thole_a is required when use_thole = true.');
        end
        thole_a = sys.thole_a;
    else
        thole_a = [];
    end

    activeSites = rowGeomCache.activeSites(:);
    srcFull = rowGeomCache.source_full_idx(:);
    row_ptr = rowGeomCache.row_ptr;
    r = rowGeomCache.r_bare;
    r2 = rowGeomCache.r2_bare;

    nInt = numel(r);
    coeff_iso = zeros(nInt, 1);
    coeff_dyad = zeros(nInt, 1);

    alpha2 = alpha^2;
    twoAlphaOverSqrtPi = 2 * alpha / sqrt(pi);

    invr2 = 1 ./ r2;
    invr = 1 ./ r;
    invr3 = invr .* invr2;
    invr5 = invr3 .* invr2;
    invr4 = invr2.^2;

    erfcar = erfc(alpha * r);
    expar2 = exp(-alpha2 * r2);

    B = erfcar .* invr3 + twoAlphaOverSqrtPi * expar2 .* invr2;
    C = 3 .* erfcar .* invr5 + twoAlphaOverSqrtPi .* (2 .* alpha2 .* invr2 + 3 .* invr4) .* expar2;

    coeff_iso = -B;
    coeff_dyad = +C;

    if use_thole
        for a = 1:numel(activeSites)
            iFull = activeSites(a);
            idx0 = row_ptr(a);
            idx1 = row_ptr(a + 1) - 1;
            if idx1 < idx0
                continue;
            end

            ai = alpha_site(iFull);
            srcIdx = srcFull(idx0:idx1);

            tf = local_thole_factors_vectorized(r(idx0:idx1), ai, alpha_site(srcIdx), thole_a);
            coeff_iso(idx0:idx1) = coeff_iso(idx0:idx1) - tf.l3 .* invr3(idx0:idx1);
            coeff_dyad(idx0:idx1) = coeff_dyad(idx0:idx1) + 3 .* tf.l5 .* invr5(idx0:idx1);
        end
    end

    coeffCache = struct();
    coeffCache.mode = 'periodic_realspace_dipole_coeff';
    coeffCache.coeff_iso = coeff_iso;
    coeffCache.coeff_dyad = coeff_dyad;
    coeffCache.alpha = alpha;
    coeffCache.rcut = rowGeomCache.rcut;
    coeffCache.use_thole = use_thole;
end

function tf = local_thole_factors_vectorized(r, alpha_i, alpha_j, thole_a)
    r = r(:);
    alpha_j = alpha_j(:);

    tf = struct();
    tf.f3 = ones(size(r));
    tf.f5 = ones(size(r));
    tf.l3 = zeros(size(r));
    tf.l5 = zeros(size(r));

    if thole_a == 0 || alpha_i == 0
        return;
    end

    keep = alpha_j > 0;
    if ~any(keep)
        return;
    end

    u = zeros(size(r));
    u(keep) = r(keep) ./ (alpha_i .* alpha_j(keep)).^(1/6);
    au3 = thole_a * u.^3;

    tf.f3(keep) = 1 - exp(-au3(keep));
    tf.f5(keep) = 1 - (1 + au3(keep)) .* exp(-au3(keep));
    tf.l3(keep) = tf.f3(keep) - 1;
    tf.l5(keep) = tf.f5(keep) - 1;
end