function [E, parts] = induced_field_from_charges_periodic(sys, fieldParams)
%INDUCED_FIELD_FROM_CHARGES_PERIODIC Periodic Ewald field from fixed charges.
%
% E = thole.induced_field_from_charges_periodic(sys, fieldParams)
% [E, parts] = thole.induced_field_from_charges_periodic(sys, fieldParams)
%
% Inputs
%   sys canonical polarization system in atomic units with fields:
%       .site_pos
%       .site_charge
%       .super_lattice or .lattice
%       .site_alpha, .thole_a   (only needed if use_thole_damping = true)
%
%   IMPORTANT LATTICE CONVENTION
%   ----------------------------
%   The builder/geom coordinate convention stores direct lattice vectors as
%   ROWS:
%
%       cart = frac * Hrow
%
%   The Ewald k-vector helper ewald.enumerate_kvecs_triclinic expects
%   direct lattice vectors as COLUMNS:
%
%       Hcol = Hrow.'
%
%   This function therefore reads the system lattice in row convention and
%   transposes it before using Ewald/periodic helper routines.
%
%   fieldParams struct with fields:
%       .exclude_self         logical, default true
%       .use_thole_damping    logical, default false
%       .target_mask          logical N x 1, default site_is_polarizable or all
%       .source_mask          logical N x 1, default abs(site_charge) > 0
%       .real_only            logical, default false
%                             if true, compute only real-space contribution
%                             and skip reciprocal/surface terms
%       .kspace_mode          'auto' | 'full' | 'blocked'
%                             legacy alias 'chunked' is accepted
%       .k_block_size         positive integer, default 2048
%       .kspace_memory_limit_gb positive scalar, default 8
%       .verbose              logical, default false
%       .ewald                struct with fields .alpha .rcut .kcut
%                             and optional .boundary ('tinfoil' or 'vacuum')
%
% Output
%   E       N x 3 field at all sites; only target rows are populated
%   parts   optional diagnostic struct with:
%             .real
%             .recip
%             .surf
%             .total
%             .Mq
%             .Esurf_q
%             .surf_coeff
%             .target_mask
%             .source_mask
%             .target_sites
%             .source_sites
%             .qtot
%             .nTargets
%             .nSources
%             .nRealEntries
%             .nK
%             .alpha
%             .rcut
%             .kcut
%             .boundary
%             .use_thole_damping
%             .exclude_self
%             .real_only
%             .storage_mode
%             .kspace_mode_requested
%             .k_block_size
%             .kspace_memory_limit_gb
%             .estimated_full_gb
%             .time_real
%             .time_recip
%             .time_total
%             .realCache
%
% Notes
% -----
% Real-space:
%   Uses geom.build_target_source_field_cache_periodic(...)
%
% Reciprocal-space:
%   Uses the standard periodic Ewald charge field. Since
%   ewald.enumerate_kvecs_triclinic returns one representative from each
%   +/- k pair, the real-valued half-k formula is
%
%     rho_cos(k) = sum_j q_j cos(k.r_j)
%     rho_sin(k) = sum_j q_j sin(k.r_j)
%
%     E_k(r_i) = sum_{k in half-space}
%                pref(k) *
%                [sin(k.r_i) rho_cos(k) - cos(k.r_i) rho_sin(k)] * k
%
% with
%
%     pref(k) = +(8*pi/V) * exp(-k^2/(4 alpha^2)) / k^2.
%
% This implementation requires net-neutral selected sources.
%
% Real-only mode:
%   fieldParams.real_only = true returns only the real-space contribution.
%   This is intended for mixed P3M workflows:
%
%       E_total = E_real,MEX/cache + E_recip,P3M + E_surface
%
%   In real-only mode, reciprocal k-vector enumeration and surface terms are
%   skipped deliberately.

    if nargin < 2 || isempty(fieldParams)
        fieldParams = struct();
    end

    wantParts = (nargout >= 2);
    tTotal = tic;

    io.assert_atomic_units(sys);

    pos = local_require_field(sys, 'site_pos');
    q   = local_require_field(sys, 'site_charge');
    q   = q(:);

    nSites = size(pos, 1);
    if numel(q) ~= nSites
        error('thole:induced_field_from_charges_periodic:BadChargeSize', ...
            'sys.site_charge must have length N.');
    end

    % System lattice is stored in row-vector convention.
    % Ewald helpers expect column-vector convention.
    Hrow = local_get_direct_lattice(sys);
    H = Hrow.';

    V = abs(det(H));
    if V <= 1e-14
        error('thole:induced_field_from_charges_periodic:SingularCell', ...
            'Direct lattice has nonpositive volume.');
    end

    exclude_self = true;
    if isfield(fieldParams, 'exclude_self') && ~isempty(fieldParams.exclude_self)
        exclude_self = logical(fieldParams.exclude_self);
    end

    useThole = false;
    if isfield(fieldParams, 'use_thole_damping') && ~isempty(fieldParams.use_thole_damping)
        useThole = logical(fieldParams.use_thole_damping);
    end

    verbose = false;
    if isfield(fieldParams, 'verbose') && ~isempty(fieldParams.verbose)
        verbose = logical(fieldParams.verbose);
    end

    real_only = false;
    if isfield(fieldParams, 'real_only') && ~isempty(fieldParams.real_only)
        real_only = logical(fieldParams.real_only);
    end

    if isfield(fieldParams, 'target_mask') && ~isempty(fieldParams.target_mask)
        target_mask = logical(fieldParams.target_mask(:));
    elseif isfield(sys, 'site_is_polarizable') && ~isempty(sys.site_is_polarizable)
        target_mask = logical(sys.site_is_polarizable(:));
    else
        target_mask = true(nSites, 1);
    end

    if isfield(fieldParams, 'source_mask') && ~isempty(fieldParams.source_mask)
        source_mask = logical(fieldParams.source_mask(:));
    else
        source_mask = abs(q) > 0;
    end

    if numel(target_mask) ~= nSites || numel(source_mask) ~= nSites
        error('thole:induced_field_from_charges_periodic:BadMask', ...
            'target_mask and source_mask must have length N.');
    end

    if ~isfield(fieldParams, 'ewald') || isempty(fieldParams.ewald)
        error('thole:induced_field_from_charges_periodic:MissingEwald', ...
            'fieldParams.ewald is required.');
    end
    ew = fieldParams.ewald;

    alpha = local_require_struct_field(ew, 'alpha');
    rcut  = local_require_struct_field(ew, 'rcut');
    kcut  = local_require_struct_field(ew, 'kcut');

    boundary = 'tinfoil';
    if isfield(ew, 'boundary') && ~isempty(ew.boundary)
        boundary = lower(char(string(ew.boundary)));
    end
    if ~ismember(boundary, {'tinfoil','vacuum'})
        error('thole:induced_field_from_charges_periodic:BoundaryNotSupported', ...
            'boundary must be ''tinfoil'' or ''vacuum''.');
    end

    kspace_mode = 'auto';
    if isfield(fieldParams, 'kspace_mode') && ~isempty(fieldParams.kspace_mode)
        kspace_mode = lower(char(string(fieldParams.kspace_mode)));
    end
    if strcmp(kspace_mode, 'chunked')
        kspace_mode = 'blocked';
    end
    if ~ismember(kspace_mode, {'auto','full','blocked'})
        error('thole:induced_field_from_charges_periodic:BadKspaceMode', ...
            'kspace_mode must be ''auto'', ''full'', or ''blocked''.');
    end

    k_block_size = 2048;
    if isfield(fieldParams, 'k_block_size') && ~isempty(fieldParams.k_block_size)
        k_block_size = fieldParams.k_block_size;
    end
    validateattributes(k_block_size, {'numeric'}, ...
        {'scalar','real','finite','positive','integer'}, ...
        mfilename, 'fieldParams.k_block_size');
    k_block_size = double(k_block_size);

    kspace_memory_limit_gb = 8;
    if isfield(fieldParams, 'kspace_memory_limit_gb') && ~isempty(fieldParams.kspace_memory_limit_gb)
        kspace_memory_limit_gb = fieldParams.kspace_memory_limit_gb;
    end
    validateattributes(kspace_memory_limit_gb, {'numeric'}, ...
        {'scalar','real','finite','positive'}, ...
        mfilename, 'fieldParams.kspace_memory_limit_gb');

    source_sites = find(source_mask);
    target_sites = find(target_mask);
    nSources = numel(source_sites);
    nTargets = numel(target_sites);

    E = zeros(nSites, 3);
    Ereal_all = zeros(nSites, 3);
    Erecip_all = zeros(nSites, 3);

    parts = struct();
    if nSources == 0 || nTargets == 0
        if wantParts
            parts = local_make_empty_parts(sys, target_mask, source_mask, ...
                target_sites, source_sites, alpha, rcut, kcut, boundary, ...
                kspace_mode, k_block_size, 0.0, toc(tTotal), real_only);
        end
        return;
    end

    qsrc = q(source_sites);
    qtot = sum(qsrc);
    if abs(qtot) > 1e-10
        error('thole:induced_field_from_charges_periodic:NonNeutralSources', ...
            ['Periodic Ewald field from fixed charges requires net-neutral selected sources. ' ...
             'Selected total charge = %+0.16e'], qtot);
    end

    if verbose
        fprintf('induced_field_from_charges_periodic:\n');
        fprintf('  nTargets = %d\n', nTargets);
        fprintf('  nSources = %d\n', nSources);
        fprintf('  alpha    = %.6f\n', alpha);
        fprintf('  rcut     = %.6f bohr\n', rcut);
        fprintf('  kcut     = %.6f bohr^-1\n', kcut);
        fprintf('  boundary = %s\n', boundary);
        fprintf('  use_thole_damping = %d\n', useThole);
        fprintf('  real_only = %d\n', real_only);
    end

    %% --------------------------------------------------------------------
    % Real-space contribution via periodic target-source cache

    tReal = tic;

    cacheOpts = struct();
    cacheOpts.exclude_self = exclude_self;
    cacheOpts.use_thole = useThole;
    cacheOpts.profile = verbose;

    realCache = geom.build_target_source_field_cache_periodic( ...
        sys, target_mask, source_mask, ew, cacheOpts);

    ErealTarget = zeros(nTargets, 3);

    row_ptr = realCache.row_ptr;
    col_idx = realCache.col_idx;
    dr      = realCache.dr;
    coeff   = realCache.inv_r3_eff;

    for it = 1:nTargets
        i0 = row_ptr(it);
        i1 = row_ptr(it + 1) - 1;

        if i1 < i0
            continue;
        end

        idx = i0:i1;
        qj = qsrc(col_idx(idx));
        c = qj .* coeff(idx);
        ErealTarget(it, :) = sum(c .* dr(idx, :), 1);
    end

    time_real = toc(tReal);

    Etarget = ErealTarget;
    Ereal_all(target_sites, :) = ErealTarget;

    %% --------------------------------------------------------------------
    % Optional real-only return.
    %
    % This is useful for mixed P3M workflows:
    %
    %   E_total = E_real,MEX/cache + E_recip,P3M + E_surf
    %
    % In real-only mode, this function deliberately does not enumerate
    % k-vectors and does not add the surface term.

    if real_only
        E(target_sites, :) = Etarget;

        time_recip = 0.0;
        nk = 0;
        storage_mode = 'real_only';
        estimated_full_gb = 0.0;

        Esurf_all = zeros(nSites, 3);
        Mq = [0.0, 0.0, 0.0];
        Esurf_q = [0.0, 0.0, 0.0];
        surf_coeff = 0.0;

        if wantParts
            parts = struct();
            parts.real = Ereal_all;
            parts.recip = Erecip_all;
            parts.surf = Esurf_all;
            parts.Mq = Mq;
            parts.Esurf_q = Esurf_q;
            parts.surf_coeff = surf_coeff;
            parts.total = E;

            parts.target_mask = target_mask;
            parts.source_mask = source_mask;
            parts.target_sites = target_sites;
            parts.source_sites = source_sites;
            parts.qtot = qtot;

            parts.nTargets = nTargets;
            parts.nSources = nSources;
            parts.nRealEntries = realCache.nEntries;
            parts.nK = nk;

            parts.alpha = alpha;
            parts.rcut = rcut;
            parts.kcut = kcut;
            parts.boundary = boundary;
            parts.use_thole_damping = useThole;
            parts.exclude_self = exclude_self;
            parts.real_only = true;

            parts.storage_mode = storage_mode;
            parts.kspace_mode_requested = kspace_mode;
            parts.k_block_size = k_block_size;
            parts.kspace_memory_limit_gb = kspace_memory_limit_gb;
            parts.estimated_full_gb = estimated_full_gb;

            parts.time_real = time_real;
            parts.time_recip = time_recip;
            parts.time_total = toc(tTotal);

            parts.realCache = realCache;
        end

        if verbose
            fprintf('  real_only mode: skipping reciprocal and surface terms\n');
            fprintf('  ||Ereal||_F  = %.16e\n', norm(Ereal_all, 'fro'));
            fprintf('  ||Etotal||_F = %.16e\n', norm(E, 'fro'));
            fprintf('  time real    = %.6f s\n', time_real);
            fprintf('  time recip   = %.6f s\n', time_recip);
            fprintf('  time total   = %.6f s\n', toc(tTotal));
        end

        return;
    end

    %% --------------------------------------------------------------------
    % Reciprocal-space contribution

    tRecip = tic;

    [kvecs, meta] = ewald.enumerate_kvecs_triclinic(H, kcut);
    nk = size(kvecs, 1);

    storage_mode = 'none';
    estimated_full_gb = 0.0;

    if nk > 0
        k2 = meta.k2(:);

        % Half-k real-valued charge-field prefactor.
        % enumerate_kvecs_triclinic keeps one representative from each +/-k
        % pair, so the full complex 4*pi/V sum becomes 8*pi/V here.
        %
        % With amp = sin(k.r_i) rho_cos - cos(k.r_i) rho_sin, the electric
        % field -grad(phi) has a positive prefactor.
        pref = +(8 * pi / V) * exp(-k2 ./ (4 * alpha^2)) ./ k2;

        source_pos = pos(source_sites, :);
        target_pos = pos(target_sites, :);

        estimated_full_bytes = double(nSources + nTargets) * double(nk) * 8 * 3;
        estimated_full_gb = estimated_full_bytes / 1024^3;
        memory_limit_bytes = kspace_memory_limit_gb * 1024^3;

        if strcmp(kspace_mode, 'auto')
            if estimated_full_bytes > memory_limit_bytes
                storage_mode = 'blocked';
            else
                storage_mode = 'full';
            end
        else
            storage_mode = kspace_mode;
        end

        if verbose
            fprintf('  reciprocal nK = %d | mode = %s | est full = %.3f GB\n', ...
                nk, storage_mode, estimated_full_gb);
        end

        ErecipTarget = zeros(nTargets, 3);

        switch storage_mode
            case 'full'
                phase_src = source_pos * kvecs.';
                cos_src = cos(phase_src);
                sin_src = sin(phase_src);

                rho_cos = qsrc.' * cos_src;   % 1 x nk
                rho_sin = qsrc.' * sin_src;   % 1 x nk

                phase_tgt = target_pos * kvecs.';
                cos_tgt = cos(phase_tgt);
                sin_tgt = sin(phase_tgt);

                amp = sin_tgt .* rho_cos - cos_tgt .* rho_sin;  % nTargets x nk
                w = amp .* pref.';                               % nTargets x nk

                ErecipTarget(:,1) = ErecipTarget(:,1) + w * kvecs(:,1);
                ErecipTarget(:,2) = ErecipTarget(:,2) + w * kvecs(:,2);
                ErecipTarget(:,3) = ErecipTarget(:,3) + w * kvecs(:,3);

            case 'blocked'
                blockStart = (1:k_block_size:nk).';
                nBlocks = numel(blockStart);

                for b = 1:nBlocks
                    i0 = blockStart(b);
                    i1 = min(i0 + k_block_size - 1, nk);
                    idx = i0:i1;

                    kb = kvecs(idx, :);
                    prefb = pref(idx);

                    phase_src = source_pos * kb.';
                    cos_src = cos(phase_src);
                    sin_src = sin(phase_src);

                    rho_cos = qsrc.' * cos_src;
                    rho_sin = qsrc.' * sin_src;

                    phase_tgt = target_pos * kb.';
                    cos_tgt = cos(phase_tgt);
                    sin_tgt = sin(phase_tgt);

                    amp = sin_tgt .* rho_cos - cos_tgt .* rho_sin;
                    w = amp .* prefb.';

                    ErecipTarget(:,1) = ErecipTarget(:,1) + w * kb(:,1);
                    ErecipTarget(:,2) = ErecipTarget(:,2) + w * kb(:,2);
                    ErecipTarget(:,3) = ErecipTarget(:,3) + w * kb(:,3);
                end

            otherwise
                error('thole:induced_field_from_charges_periodic:InternalModeError', ...
                    'Unexpected storage mode.');
        end

        Etarget = Etarget + ErecipTarget;
        Erecip_all(target_sites, :) = ErecipTarget;
    end

    time_recip = toc(tRecip);

    %% --------------------------------------------------------------------
    % Surface term from fixed-charge cell dipole

    Esurf_all = zeros(nSites, 3);
    Mq = [0.0, 0.0, 0.0];
    Esurf_q = [0.0, 0.0, 0.0];
    surf_coeff = 0.0;

    switch boundary
        case 'tinfoil'
            % no surface term

        case 'vacuum'
            % Neutral source charge distribution, so Mq is origin-independent
            % for a fixed Cartesian branch.
            Mq = sum(qsrc .* pos(source_sites, :), 1);

            surf_coeff = 4 * pi / (3 * V);
            Esurf_q = -surf_coeff * Mq;

            Etarget = Etarget + Esurf_q;
            Esurf_all(target_sites, :) = repmat(Esurf_q, nTargets, 1);

        otherwise
            error('thole:induced_field_from_charges_periodic:UnknownBoundary', ...
                'Unknown boundary "%s".', boundary);
    end

    E(target_sites, :) = Etarget;

    if wantParts
        parts = struct();
        parts.real = Ereal_all;
        parts.recip = Erecip_all;
        parts.surf = Esurf_all;
        parts.Mq = Mq;
        parts.Esurf_q = Esurf_q;
        parts.surf_coeff = surf_coeff;
        parts.total = E;

        parts.target_mask = target_mask;
        parts.source_mask = source_mask;
        parts.target_sites = target_sites;
        parts.source_sites = source_sites;
        parts.qtot = qtot;

        parts.nTargets = nTargets;
        parts.nSources = nSources;
        parts.nRealEntries = realCache.nEntries;
        parts.nK = nk;

        parts.alpha = alpha;
        parts.rcut = rcut;
        parts.kcut = kcut;
        parts.boundary = boundary;
        parts.use_thole_damping = useThole;
        parts.exclude_self = exclude_self;
        parts.real_only = false;

        parts.storage_mode = storage_mode;
        parts.kspace_mode_requested = kspace_mode;
        parts.k_block_size = k_block_size;
        parts.kspace_memory_limit_gb = kspace_memory_limit_gb;
        parts.estimated_full_gb = estimated_full_gb;

        parts.time_real = time_real;
        parts.time_recip = time_recip;
        parts.time_total = toc(tTotal);

        parts.realCache = realCache;
    end

    if verbose
        fprintf('  ||Ereal||_F  = %.16e\n', norm(Ereal_all, 'fro'));
        fprintf('  ||Erecip||_F = %.16e\n', norm(Erecip_all, 'fro'));
        fprintf('  ||Esurf||_F  = %.16e\n', norm(Esurf_all, 'fro'));
        fprintf('  ||Etotal||_F = %.16e\n', norm(E, 'fro'));
        fprintf('  time real    = %.6f s\n', time_real);
        fprintf('  time recip   = %.6f s\n', time_recip);
        fprintf('  time total   = %.6f s\n', toc(tTotal));
    end
end

function parts = local_make_empty_parts(sys, target_mask, source_mask, ...
    target_sites, source_sites, alpha, rcut, kcut, boundary, ...
    kspace_mode, k_block_size, time_before_return, time_total, real_only)

    nSites = size(sys.site_pos, 1);

    parts = struct();
    parts.real = zeros(nSites, 3);
    parts.recip = zeros(nSites, 3);
    parts.surf = zeros(nSites, 3);
    parts.Mq = [0.0, 0.0, 0.0];
    parts.Esurf_q = [0.0, 0.0, 0.0];
    parts.surf_coeff = 0.0;
    parts.total = zeros(nSites, 3);

    parts.target_mask = target_mask;
    parts.source_mask = source_mask;
    parts.target_sites = target_sites;
    parts.source_sites = source_sites;
    parts.qtot = 0.0;

    parts.nTargets = numel(target_sites);
    parts.nSources = numel(source_sites);
    parts.nRealEntries = 0;
    parts.nK = 0;

    parts.alpha = alpha;
    parts.rcut = rcut;
    parts.kcut = kcut;
    parts.boundary = boundary;
    parts.use_thole_damping = false;
    parts.exclude_self = true;
    parts.real_only = real_only;

    if real_only
        parts.storage_mode = 'real_only';
    else
        parts.storage_mode = 'none';
    end
    parts.kspace_mode_requested = kspace_mode;
    parts.k_block_size = k_block_size;
    parts.kspace_memory_limit_gb = NaN;
    parts.estimated_full_gb = 0.0;

    parts.time_real = time_before_return;
    parts.time_recip = 0.0;
    parts.time_total = time_total;

    parts.realCache = struct();
end

function x = local_require_field(s, name)
    if ~isfield(s, name) || isempty(s.(name))
        error('thole:induced_field_from_charges_periodic:MissingField', ...
            'sys.%s is required.', name);
    end
    x = s.(name);
end

function x = local_require_struct_field(s, name)
    if ~isfield(s, name) || isempty(s.(name))
        error('thole:induced_field_from_charges_periodic:MissingEwaldField', ...
            'fieldParams.ewald.%s is required.', name);
    end
    x = s.(name);
end

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('thole:induced_field_from_charges_periodic:MissingLattice', ...
            'sys.super_lattice or sys.lattice is required.');
    end

    validateattributes(H, {'double'}, {'size',[3,3],'real','finite'}, ...
        mfilename, 'lattice');
end