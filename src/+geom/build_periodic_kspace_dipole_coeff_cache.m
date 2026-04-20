function coeffCache = build_periodic_kspace_dipole_coeff_cache(targetCache)
%BUILD_PERIODIC_KSPACE_DIPOLE_COEFF_CACHE Build dipole-specific reciprocal payload.
%
% coeffCache = geom.build_periodic_kspace_dipole_coeff_cache(targetCache)

    narginchk(1, 1);

    nk = targetCache.num_kvec;
    kvecs = targetCache.kvecs;

    kx = kvecs(:, 1);
    ky = kvecs(:, 2);
    kz = kvecs(:, 3);

    kk6 = zeros(nk, 6);
    kk6(:, 1) = kx .* kx;
    kk6(:, 2) = ky .* ky;
    kk6(:, 3) = kz .* kz;
    kk6(:, 4) = kx .* ky;
    kk6(:, 5) = kx .* kz;
    kk6(:, 6) = ky .* kz;

    % Match old build_periodic_kspace_cache exactly
    pref = -targetCache.pref_base;

    activeSites = targetCache.targetSites(:);
    nSites = targetCache.nSites;
    fullToActive = zeros(nSites, 1);
    fullToActive(activeSites) = 1:numel(activeSites);

    coeffCache = struct();
    coeffCache.mode = 'periodic_kspace_dipole_coeff';
    coeffCache.pref = pref;
    coeffCache.kk6 = kk6;
    coeffCache.activeSites = activeSites;
    coeffCache.full_to_active = fullToActive;
end