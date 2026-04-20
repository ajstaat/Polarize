function coeffCache = build_periodic_kspace_charge_coeff_cache(targetCache, sourceCache, sourceCharge)
%BUILD_PERIODIC_KSPACE_CHARGE_COEFF_CACHE Build charge-specific reciprocal payload.
%
% coeffCache = geom.build_periodic_kspace_charge_coeff_cache(targetCache, sourceCache, sourceCharge)

    narginchk(3, 3);

    sourceCharge = sourceCharge(:);
    if numel(sourceCharge) ~= numel(sourceCache.sourceSites)
        error('geom:build_periodic_kspace_charge_coeff_cache:BadChargeSize', ...
            'sourceCharge must have length equal to numel(sourceCache.sourceSites).');
    end

    Cq = (sourceCharge.' * sourceCache.source_cos).';   % Nk x 1
    Sq = (sourceCharge.' * sourceCache.source_sin).';   % Nk x 1

    coeffCache = struct();
    coeffCache.mode = 'periodic_kspace_charge_coeff';
    coeffCache.pref = 2 * targetCache.pref_base;   % half-space enumeration
    coeffCache.sourceSites = sourceCache.sourceSites(:);
    coeffCache.source_charge = sourceCharge;
    coeffCache.rho_cos = Cq;
    coeffCache.rho_sin = Sq;
end