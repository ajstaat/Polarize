function Edip = induced_field_from_dipoles_thole_periodic(sys, mu, ewaldParams, dipoleParams)
%INDUCED_FIELD_FROM_DIPOLES_THOLE_PERIODIC Periodic induced field from dipoles.
%
% Edip = thole.induced_field_from_dipoles_thole_periodic(sys, mu, ewaldParams)
% Edip = thole.induced_field_from_dipoles_thole_periodic(sys, mu, ewaldParams, dipoleParams)
%
% Inputs
%   sys canonical polarization-system struct in atomic units
%   mu  N x 3 induced dipoles
%   ewaldParams struct with fields:
%       .alpha
%       .rcut
%       .kcut
%       .boundary   optional, default 'tinfoil'
%
%   dipoleParams optional struct with fields:
%       .target_mask             N x 1 logical, optional
%       .source_mask             N x 1 logical, optional
%       .use_thole               logical, default true
%       .realspace_row_cache     optional active periodic row cache from
%                                geom.build_active_row_cache_periodic(...)
%       .kspace_cache            optional cache from
%                                geom.build_periodic_kspace_cache(...)
%       .problem                 optional problem struct from
%                                thole.prepare_scf_problem(...)
%       .kspace_mode             'auto' | 'full' | 'blocked'
%                                (legacy alias 'chunked' accepted)
%       .kspace_memory_limit_gb  positive scalar, default 8
%       .k_block_size            positive integer, default 2048
%       .verbose                 logical, default false
%
% Output
%   Edip N x 3 periodic induced field
%
% Notes
% -----
% - This version assumes the minimal-branch cache set:
%     * real-space uses only the active periodic row cache
%     * reciprocal-space uses full or blocked k-space cache
% - Legacy unordered / layered periodic real-space cache paths have been removed.
% - The active periodic row cache is interpreted in the active-space basis:
%     row_ptr / col_idx / dr / coeff_iso / coeff_dyad
%   with col_idx indexing activeSites.

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
    end

    % ---------------------------------------------------------------------
    % Real-space active periodic row cache

    if isfield(dipoleParams, 'realspace_row_cache') && ~isempty(dipoleParams.realspace_row_cache)
        rowCache = dipoleParams.realspace_row_cache;
    else
        rowOpts = struct();
        rowOpts.profile   = verbose;
        rowOpts.use_mex   = false;
        rowOpts.use_thole = use_thole;
        rowCache = geom.build_active_row_cache_periodic(sys, problem, ewaldParams, rowOpts);
    end

    % ---------------------------------------------------------------------
    % K-space cache

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

    if strcmpi(kCache.storage_mode, 'chunked')
        kCache.storage_mode = 'blocked';
    end
    if isfield(kCache, 'phase_storage') && strcmpi(kCache.phase_storage, 'chunked')
        kCache.phase_storage = 'blocked';
    end

    Edip = zeros(nSites, 3);

    % ---------------------------------------------------------------------
    % Real-space contribution from active periodic row cache

    row_ptr = rowCache.row_ptr(:);
    col_idx = rowCache.col_idx(:);
    dr_all  = rowCache.dr;

    if ~isfield(rowCache, 'coeff_iso') || ~isfield(rowCache, 'coeff_dyad')
        error('thole:induced_field_from_dipoles_thole_periodic:MissingRowCoeffs', ...
            'realspace_row_cache must provide coeff_iso and coeff_dyad.');
    end
    coeff_iso  = rowCache.coeff_iso(:);
    coeff_dyad = rowCache.coeff_dyad(:);

    if isfield(rowCache, 'activeSites') && ~isempty(rowCache.activeSites)
        activeSites = rowCache.activeSites(:);
    elseif isfield(problem, 'activeSites') && ~isempty(problem.activeSites)
        activeSites = problem.activeSites(:);
    else
        error('thole:induced_field_from_dipoles_thole_periodic:MissingActiveSites', ...
            'Need activeSites from row cache or problem.');
    end

    nActive = numel(activeSites);

    for a = 1:nActive
        iFull = activeSites(a);
        if ~target_mask(iFull)
            continue;
        end

        idx0 = row_ptr(a);
        idx1 = row_ptr(a + 1) - 1;
        if idx1 < idx0
            continue;
        end

        idx = idx0:idx1;
        srcActive = col_idx(idx);
        srcFull = activeSites(srcActive);

        keep = source_mask(srcFull);
        if ~any(keep)
            continue;
        end

        idx = idx(keep);
        srcFull = srcFull(keep);

        muSrc = mu(srcFull, :);
        dr = dr_all(idx, :);

        muDotR = sum(muSrc .* dr, 2);
        contrib = coeff_iso(idx) .* muSrc + coeff_dyad(idx) .* (muDotR .* dr);

        Edip(iFull, :) = Edip(iFull, :) + sum(contrib, 1);
    end

    % ---------------------------------------------------------------------
    % Reciprocal-space contribution

    activeSites = kCache.activeSites(:);
    nPol = kCache.nPolSites;

    if kCache.num_kvec > 0 && nPol > 0
        mu_pol = mu(activeSites, :);

        source_active = source_mask(activeSites);
        target_active = target_mask(activeSites);

        mu_src = mu_pol;
        mu_src(~source_active, :) = 0;

        switch lower(kCache.storage_mode)
            case 'full'
                kvecs = kCache.kvecs;
                two_pref = kCache.two_pref(:);

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

            case 'blocked'
                pos_pol = kCache.active_pos;
                Erecip_pol = zeros(nPol, 3);

                for b = 1:kCache.num_blocks
                    blk = kCache.blocks(b);
                    if blk.nk == 0
                        continue;
                    end

                    phase_blk = pos_pol * blk.kvecs.';
                    cos_blk = cos(phase_blk);
                    sin_blk = sin(phase_blk);

                    v_blk = mu_src * blk.kvecs.';
                    A_blk = sum(cos_blk .* v_blk, 1);
                    B_blk = sum(sin_blk .* v_blk, 1);

                    phase_factor_blk = cos_blk .* A_blk + sin_blk .* B_blk;
                    W_blk = phase_factor_blk .* (blk.two_pref.');

                    Erecip_pol(:, 1) = Erecip_pol(:, 1) + W_blk * blk.kvecs(:, 1);
                    Erecip_pol(:, 2) = Erecip_pol(:, 2) + W_blk * blk.kvecs(:, 2);
                    Erecip_pol(:, 3) = Erecip_pol(:, 3) + W_blk * blk.kvecs(:, 3);
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
            surf_coeff = -4 * pi / (3 * V);

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