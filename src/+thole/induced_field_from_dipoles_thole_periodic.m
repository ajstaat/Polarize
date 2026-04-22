function Edip = induced_field_from_dipoles_thole_periodic(sys, mu, ewaldParams, dipoleParams)
%INDUCED_FIELD_FROM_DIPOLES_THOLE_PERIODIC Periodic induced field from dipoles.
%
% Edip = thole.induced_field_from_dipoles_thole_periodic(sys, mu, ewaldParams)
% Edip = thole.induced_field_from_dipoles_thole_periodic(sys, mu, ewaldParams, dipoleParams)
%
% Inputs
%   sys           canonical polarization-system struct in atomic units
%   mu            N x 3 induced dipoles
%   ewaldParams   struct with fields:
%                   .alpha
%                   .rcut
%                   .kcut
%                   .boundary   optional, default 'tinfoil'
%
%   dipoleParams  optional struct with fields:
%                   .target_mask                 N x 1 logical, optional
%                   .source_mask                 N x 1 logical, optional
%                   .use_thole                   logical, default true
%                   .realspace_cache             optional LEGACY unordered cache
%                   .realspace_row_cache         optional directed row cache
%                   .realspace_dipole_coeff_cache optional coeff cache for
%                                                layered row-geometry path
%                   .kspace_cache                optional cache from
%                                                geom.build_periodic_kspace_cache(...)
%                   .problem                     optional problem struct from
%                                                thole.prepare_scf_problem(...)
%                   .kspace_mode                 'auto' | 'full' | 'chunked'
%                   .kspace_memory_limit_gb      positive scalar, default 8
%                   .k_block_size                positive integer, default 2048
%                   .verbose                     logical, default false
%
% Output
%   Edip          N x 3 periodic induced field
%
% Real-space supported modes
%   (A) Legacy unordered cache:
%         pair_i, pair_j, dr, coeff_iso, coeff_dyad, nInteractions
%
%   (B) Combined directed row cache:
%         row_ptr, source_full_idx, dr, coeff_iso, coeff_dyad
%         plus targetSites or activeSites
%
%   (C) Layered row cache:
%         realspace_row_cache            (geometry)
%         realspace_dipole_coeff_cache   (coeff_iso / coeff_dyad)

    io.assert_atomic_units(sys);

    if nargin < 4 || isempty(dipoleParams)
        dipoleParams = struct();
    end

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('thole:induced_field_from_dipoles_thole_periodic:MissingSitePos', ...
            'sys.site_pos is missing or empty.');
    end
    if ~isfield(sys, 'site_is_polarizable') || isempty(sys.site_is_polarizable)
        error('thole:induced_field_from_dipoles_thole_periodic:MissingSiteMask', ...
            'sys.site_is_polarizable is missing or empty.');
    end
    if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
        error('thole:induced_field_from_dipoles_thole_periodic:MissingSiteAlpha', ...
            'sys.site_alpha is missing or empty.');
    end

    nSites = size(sys.site_pos, 1);
    if ~isequal(size(mu), [nSites, 3])
        error('thole:induced_field_from_dipoles_thole_periodic:BadMuSize', ...
            'mu must be N x 3.');
    end

    if ~isfield(ewaldParams, 'alpha') || isempty(ewaldParams.alpha)
        error('thole:induced_field_from_dipoles_thole_periodic:MissingAlpha', ...
            'ewaldParams.alpha is required.');
    end
    if ~isfield(ewaldParams, 'rcut') || isempty(ewaldParams.rcut)
        error('thole:induced_field_from_dipoles_thole_periodic:MissingRcut', ...
            'ewaldParams.rcut is required.');
    end
    if ~isfield(ewaldParams, 'kcut') || isempty(ewaldParams.kcut)
        error('thole:induced_field_from_dipoles_thole_periodic:MissingKcut', ...
            'ewaldParams.kcut is required.');
    end

    target_mask = true(nSites, 1);
    if isfield(dipoleParams, 'target_mask') && ~isempty(dipoleParams.target_mask)
        target_mask = logical(dipoleParams.target_mask(:));
    end
    if numel(target_mask) ~= nSites
        error('thole:induced_field_from_dipoles_thole_periodic:BadTargetMask', ...
            'target_mask must have length N.');
    end

    source_mask = logical(sys.site_is_polarizable(:));
    if isfield(dipoleParams, 'source_mask') && ~isempty(dipoleParams.source_mask)
        source_mask = logical(dipoleParams.source_mask(:));
    end
    if numel(source_mask) ~= nSites
        error('thole:induced_field_from_dipoles_thole_periodic:BadSourceMask', ...
            'source_mask must have length N.');
    end

    use_thole = true;
    if isfield(dipoleParams, 'use_thole') && ~isempty(dipoleParams.use_thole)
        use_thole = logical(dipoleParams.use_thole);
    end

    verbose = false;
    if isfield(dipoleParams, 'verbose') && ~isempty(dipoleParams.verbose)
        verbose = logical(dipoleParams.verbose);
    end

    problem = [];
    if isfield(dipoleParams, 'problem') && ~isempty(dipoleParams.problem)
        problem = dipoleParams.problem;
    end

    if isempty(problem)
        polMask = logical(sys.site_is_polarizable(:));
        activeSites = find(polMask);
        problem = struct();
        problem.nSites = nSites;
        problem.polMask = polMask;
        problem.activeSites = activeSites;
        problem.nPolSites = numel(activeSites);
        problem.activeVecIdx = zeros(3 * numel(activeSites), 1); %#ok<STRNU>
    end

    % ---------------------------------------------------------------------
    % Real-space cache selection

    rowCache = [];
    coeffCache = [];
    haveDirectedRowCache = false;

    if isfield(dipoleParams, 'realspace_row_cache') && ~isempty(dipoleParams.realspace_row_cache)
        rowCache = dipoleParams.realspace_row_cache;
        haveDirectedRowCache = true;

        if isfield(dipoleParams, 'realspace_dipole_coeff_cache') && ...
           ~isempty(dipoleParams.realspace_dipole_coeff_cache)
            coeffCache = dipoleParams.realspace_dipole_coeff_cache;
        end
    end

    realCache = [];
    if ~haveDirectedRowCache
        if isfield(dipoleParams, 'realspace_cache') && ~isempty(dipoleParams.realspace_cache)
            realCache = dipoleParams.realspace_cache;
        else
            cacheParams = struct();
            cacheParams.use_thole = use_thole;
            cacheParams.verbose = verbose;
            realCache = geom.build_periodic_realspace_cache(sys, problem, ewaldParams, cacheParams);
        end
    end

    % ---------------------------------------------------------------------
    % K-space cache

    kCache = [];
    if isfield(dipoleParams, 'kspace_cache') && ~isempty(dipoleParams.kspace_cache)
        kCache = dipoleParams.kspace_cache;
    else
        kOpts = struct();
        if isfield(dipoleParams, 'kspace_mode') && ~isempty(dipoleParams.kspace_mode)
            kOpts.kspace_mode = dipoleParams.kspace_mode;
        end
        if isfield(dipoleParams, 'kspace_memory_limit_gb') && ~isempty(dipoleParams.kspace_memory_limit_gb)
            kOpts.kspace_memory_limit_gb = dipoleParams.kspace_memory_limit_gb;
        end
        if isfield(dipoleParams, 'k_block_size') && ~isempty(dipoleParams.k_block_size)
            kOpts.k_block_size = dipoleParams.k_block_size;
        end
        kOpts.verbose = verbose;
        kCache = geom.build_periodic_kspace_cache(sys, problem, ewaldParams, kOpts);
    end

    Edip = zeros(nSites, 3);

    % ---------------------------------------------------------------------
    % Real-space contribution

    if haveDirectedRowCache
        row_ptr = rowCache.row_ptr(:);
        srcFull = rowCache.source_full_idx(:);
        dr_all = rowCache.dr;

        if isfield(rowCache, 'coeff_iso') && isfield(rowCache, 'coeff_dyad')
            coeff_iso = rowCache.coeff_iso(:);
            coeff_dyad = rowCache.coeff_dyad(:);
        else
            if isempty(coeffCache)
                error('thole:induced_field_from_dipoles_thole_periodic:MissingRowCoeffs', ...
                    ['Directed realspace_row_cache provided without coeff_iso/coeff_dyad, ' ...
                     'and no realspace_dipole_coeff_cache was supplied.']);
            end
            coeff_iso = coeffCache.coeff_iso(:);
            coeff_dyad = coeffCache.coeff_dyad(:);
        end

        if isfield(rowCache, 'targetSites') && ~isempty(rowCache.targetSites)
            targetSites = rowCache.targetSites(:);
        elseif isfield(rowCache, 'activeSites') && ~isempty(rowCache.activeSites)
            targetSites = rowCache.activeSites(:);
        else
            error('thole:induced_field_from_dipoles_thole_periodic:MissingRowTargets', ...
                'Directed row cache must provide targetSites or activeSites.');
        end

        nTarget = numel(targetSites);

        for a = 1:nTarget
            iFull = targetSites(a);
            if ~target_mask(iFull)
                continue;
            end

            idx0 = row_ptr(a);
            idx1 = row_ptr(a + 1) - 1;
            if idx1 < idx0
                continue;
            end

            idx = idx0:idx1;
            src = srcFull(idx);

            keep = source_mask(src);
            if ~any(keep)
                continue;
            end

            idx = idx(keep);
            src = src(keep);

            muSrc = mu(src, :);
            dr = dr_all(idx, :);

            muDotR = sum(muSrc .* dr, 2);
            contrib = coeff_iso(idx) .* muSrc + ...
                      coeff_dyad(idx) .* (muDotR .* dr);

            Edip(iFull, :) = Edip(iFull, :) + sum(contrib, 1);
        end

    else
        % Legacy unordered cache path
        pair_i = realCache.pair_i(:);
        pair_j = realCache.pair_j(:);
        dr = realCache.dr;
        coeff_iso = realCache.coeff_iso(:);
        coeff_dyad = realCache.coeff_dyad(:);

        nReal = realCache.nInteractions;

        for p = 1:nReal
            i = pair_i(p);
            j = pair_j(p);

            if i == j
                if ~(target_mask(i) && source_mask(i))
                    continue;
                end

                muj = mu(i, :);
                if all(muj == 0)
                    continue;
                end

                rij = dr(p, :);
                muDotR = dot(muj, rij);

                Edip(i, :) = Edip(i, :) ...
                    + coeff_iso(p) * muj ...
                    + coeff_dyad(p) * rij * muDotR;
            else
                rij = dr(p, :);

                if target_mask(i) && source_mask(j)
                    muj = mu(j, :);
                    if ~all(muj == 0)
                        muDotR = dot(muj, rij);
                        Edip(i, :) = Edip(i, :) ...
                            + coeff_iso(p) * muj ...
                            + coeff_dyad(p) * rij * muDotR;
                    end
                end

                if target_mask(j) && source_mask(i)
                    mui = mu(i, :);
                    if ~all(mui == 0)
                        rji = -rij;
                        muDotR = dot(mui, rji);
                        Edip(j, :) = Edip(j, :) ...
                            + coeff_iso(p) * mui ...
                            + coeff_dyad(p) * rji * muDotR;
                    end
                end
            end
        end
    end

    % ---------------------------------------------------------------------
    % Reciprocal-space contribution

    activeSites = kCache.activeSites(:);
    nPol = kCache.nPolSites;

    if kCache.num_kvec > 0 && nPol > 0
        mu_pol = mu(activeSites, :);

        source_active = source_mask(activeSites);
        mu_src = mu_pol;
        mu_src(~source_active, :) = 0;

        target_active = target_mask(activeSites);

        kvecs = kCache.kvecs;
        pref = kCache.pref(:);
        nk = kCache.num_kvec;
        two_pref = 2 * pref;

        switch kCache.storage_mode
            case 'full'
                cos_phase = kCache.cos_phase;
                sin_phase = kCache.sin_phase;

                v = mu_src * kvecs.';
                A = sum(cos_phase .* v, 1);
                B = sum(sin_phase .* v, 1);

                phase_factor = cos_phase .* A + sin_phase .* B;
                W = phase_factor .* (two_pref.');

                Ex_pol = W * kvecs(:, 1);
                Ey_pol = W * kvecs(:, 2);
                Ez_pol = W * kvecs(:, 3);

                Erecip_pol = [Ex_pol, Ey_pol, Ez_pol];
                Erecip_pol(~target_active, :) = 0;

            case 'chunked'
                pos_pol = kCache.active_pos;
                blk = kCache.k_block_size;

                Erecip_pol = zeros(nPol, 3);

                for k0 = 1:blk:nk
                    k1 = min(k0 + blk - 1, nk);
                    idx = k0:k1;

                    kblk = kvecs(idx, :);
                    phase_blk = pos_pol * kblk.';
                    cos_blk = cos(phase_blk);
                    sin_blk = sin(phase_blk);

                    v_blk = mu_src * kblk.';

                    A_blk = sum(cos_blk .* v_blk, 1);
                    B_blk = sum(sin_blk .* v_blk, 1);

                    phase_factor_blk = cos_blk .* A_blk + sin_blk .* B_blk;
                    W_blk = phase_factor_blk .* (two_pref(idx).');

                    Erecip_pol(:, 1) = Erecip_pol(:, 1) + W_blk * kblk(:, 1);
                    Erecip_pol(:, 2) = Erecip_pol(:, 2) + W_blk * kblk(:, 2);
                    Erecip_pol(:, 3) = Erecip_pol(:, 3) + W_blk * kblk(:, 3);
                end

                Erecip_pol(~target_active, :) = 0;

            otherwise
                error('thole:induced_field_from_dipoles_thole_periodic:BadKspaceMode', ...
                    'Unknown kCache.storage_mode "%s".', kCache.storage_mode);
        end

        Edip(activeSites, :) = Edip(activeSites, :) + Erecip_pol;
    end

    % ---------------------------------------------------------------------
    % Analytic self term

    alpha = ewaldParams.alpha;
    self_coeff = +(4 * alpha^3 / (3 * sqrt(pi)));

    self_sites = target_mask & source_mask;
    Edip(self_sites, :) = Edip(self_sites, :) + self_coeff * mu(self_sites, :);

    % ---------------------------------------------------------------------
    % Surface term

    boundary = 'tinfoil';
    if isfield(ewaldParams, 'boundary') && ~isempty(ewaldParams.boundary)
        boundary = lower(char(string(ewaldParams.boundary)));
    end

    switch boundary
        case 'tinfoil'
            % no-op
        case 'vacuum'
            H = local_get_direct_lattice(sys);
            V = abs(det(H));
            surf_coeff = 4 * pi / (3 * V);

            Msrc = sum(mu(source_mask, :), 1);
            Edip(target_mask, :) = Edip(target_mask, :) + surf_coeff * Msrc;
        otherwise
            error('thole:induced_field_from_dipoles_thole_periodic:UnknownBoundary', ...
                'boundary must be ''tinfoil'' or ''vacuum''.');
    end
end

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('thole:induced_field_from_dipoles_thole_periodic:MissingLattice', ...
            'Need sys.super_lattice or sys.lattice.');
    end

    validateattributes(H, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
        mfilename, 'lattice');
end