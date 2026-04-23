function cache = build_target_source_field_cache_periodic(sys, target_mask, source_mask, ewaldParams, opts)
%BUILD_TARGET_SOURCE_FIELD_CACHE_PERIODIC Periodic target-source real-space cache.
%
% cache = geom.build_target_source_field_cache_periodic(sys, target_mask, source_mask, ewaldParams)
% cache = geom.build_target_source_field_cache_periodic(..., opts)
%
% Inputs
%   sys          canonical polarization-system struct in atomic units with:
%                  .site_pos       N x 3 Cartesian positions (bohr)
%                  .super_lattice  3 x 3 direct lattice, ROWS are lattice vectors
%                                  (or .lattice if .super_lattice absent)
%                  .site_alpha     N x 1 polarizabilities, needed only if
%                                  opts.use_thole = true
%                  .thole_a        scalar, needed only if opts.use_thole = true
%
%   target_mask  logical N x 1
%   source_mask  logical N x 1
%   ewaldParams  struct with:
%                  .alpha
%                  .rcut
%
%   opts         optional struct with fields:
%                  .exclude_self       logical, default true
%                  .use_thole          logical, default false
%                  .profile            logical, default false
%                  .target_block_size  positive integer, default 512
%
% Output
%   cache struct with fields:
%     .mode
%     .nTargets
%     .nSources
%     .target_sites
%     .source_sites
%     .target_pos
%     .source_pos
%     .target_alpha
%     .source_alpha
%     .row_ptr         (nTargets+1) x 1 CSR pointer
%     .col_idx         nEntries x 1 source-local column indices
%     .dr              nEntries x 3 displacement from source image to target
%     .r_bare          nEntries x 1
%     .r2_bare         nEntries x 1
%     .inv_r3_eff      nEntries x 1 effective scalar field coefficient so that
%                      E = sum_j q_j * inv_r3_eff .* dr
%     .B_ewald         nEntries x 1 Ewald real-space charge-field coefficient
%     .thole_delta     nEntries x 1 short-range Thole correction
%     .nEntries
%     .alpha
%     .rcut
%     .Lmin
%     .exclude_self
%     .use_thole
%
% Notes
% -----
% - This is an asymmetric target-source cache.  It does NOT use
%   geom.query_pairs_within_cutoff because that routine currently only
%   supports source-source/subset unordered pair queries, not arbitrary
%   target_points.
%
% - The builder/geom coordinate convention stores direct lattice vectors as
%   ROWS:
%
%       cart = frac * Hrow
%
%   For minimum-image algebra here, we convert to COLUMN-vector convention:
%
%       Hcol = Hrow.'
%
%   and use
%
%       frac = inv(Hcol) * dr_cart
%       dr_min = Hcol * (frac - round(frac)).
%
% - This cache assumes the same single-image regime as the active periodic
%   row cache:
%
%       rcut < Lmin/2
%
%   where Lmin is the shortest nonzero direct-lattice translation.
%
% - The real-space periodic Ewald charge-field contribution is
%
%       E_i^real = sum_j q_j * B_alpha(r_ij) * r_ij
%
%   with r_ij = r_i - r_j_image and
%
%       B_alpha(r) = erfc(alpha r)/r^3
%                  + (2 alpha/sqrt(pi)) exp(-alpha^2 r^2)/r^2.
%
% - If opts.use_thole = true, a short-range Thole correction is added:
%
%       B_eff(r) = B_alpha(r) + (f_Thole(r)-1)/r^3
%
%   where
%
%       f_Thole(r) = 1 - exp(-a u^3)
%       u = r / (alpha_i alpha_j)^(1/6).

    if nargin < 5 || isempty(opts)
        opts = struct();
    end

    io.assert_atomic_units(sys);

    pos = local_require_field(sys, 'site_pos');

    % Stored system lattice is row-vector convention.
    Hrow = local_get_direct_lattice(sys);

    % Internal minimum-image convention is column-vector lattice.
    Hcol = Hrow.';
    invHcol = inv(Hcol);

    nSites = size(pos, 1);

    target_mask = logical(target_mask(:));
    source_mask = logical(source_mask(:));

    if numel(target_mask) ~= nSites || numel(source_mask) ~= nSites
        error('geom:build_target_source_field_cache_periodic:BadMaskSize', ...
            'target_mask and source_mask must have length N.');
    end

    alphaEwald = local_require_struct_field(ewaldParams, 'alpha');
    rcut       = local_require_struct_field(ewaldParams, 'rcut');

    validateattributes(alphaEwald, {'double'}, {'scalar','real','positive','finite'}, ...
        mfilename, 'ewaldParams.alpha');

    validateattributes(rcut, {'double'}, {'scalar','real','positive','finite'}, ...
        mfilename, 'ewaldParams.rcut');

    exclude_self = true;
    if isfield(opts, 'exclude_self') && ~isempty(opts.exclude_self)
        exclude_self = logical(opts.exclude_self);
    end

    use_thole = false;
    if isfield(opts, 'use_thole') && ~isempty(opts.use_thole)
        use_thole = logical(opts.use_thole);
    end

    profile = false;
    if isfield(opts, 'profile') && ~isempty(opts.profile)
        profile = logical(opts.profile);
    end

    target_block_size = 512;
    if isfield(opts, 'target_block_size') && ~isempty(opts.target_block_size)
        target_block_size = opts.target_block_size;
    end
    validateattributes(target_block_size, {'numeric'}, ...
        {'scalar','integer','positive','finite'}, ...
        mfilename, 'opts.target_block_size');
    target_block_size = double(target_block_size);

    target_sites = find(target_mask);
    source_sites = find(source_mask);

    nTargets = numel(target_sites);
    nSources = numel(source_sites);

    target_pos = pos(target_sites, :);
    source_pos = pos(source_sites, :);

    if use_thole
        alpha_all = local_require_field(sys, 'site_alpha');
        alpha_all = alpha_all(:);

        if numel(alpha_all) ~= nSites
            error('geom:build_target_source_field_cache_periodic:BadAlphaSize', ...
                'sys.site_alpha must have length N.');
        end

        if ~isfield(sys, 'thole_a') || isempty(sys.thole_a)
            error('geom:build_target_source_field_cache_periodic:MissingTholeA', ...
                'sys.thole_a is required when opts.use_thole = true.');
        end

        target_alpha = alpha_all(target_sites);
        source_alpha = alpha_all(source_sites);
        thole_a = sys.thole_a;
    else
        target_alpha = zeros(nTargets, 1);
        source_alpha = zeros(nSources, 1);
        thole_a = [];
    end

    Lmin = geom.shortest_lattice_translation(Hcol);

    tol = 1e-12 * max(1, Lmin);
    rcutMaxSafe = 0.5 * Lmin - tol;

    if ~(rcut < rcutMaxSafe)
        error('geom:build_target_source_field_cache_periodic:RcutTooLarge', ...
            ['Periodic target-source field cache assumes single-image treatment.\n' ...
             'Require rcut < Lmin/2.\n' ...
             '  Lmin          = %.16g\n' ...
             '  Lmin/2        = %.16g\n' ...
             '  max safe rcut = %.16g\n' ...
             '  current rcut  = %.16g'], ...
             Lmin, 0.5 * Lmin, rcutMaxSafe, rcut);
    end

    if profile
        tStart = tic;
        fprintf('build_target_source_field_cache_periodic:\n');
        fprintf('  nTargets = %d\n', nTargets);
        fprintf('  nSources = %d\n', nSources);
        fprintf('  alpha    = %.6f bohr^-1\n', alphaEwald);
        fprintf('  rcut     = %.6f bohr\n', rcut);
        fprintf('  Lmin     = %.16e bohr\n', Lmin);
        fprintf('  use_thole = %d\n', use_thole);
        fprintf('  target_block_size = %d\n', target_block_size);
    end

    if nTargets == 0 || nSources == 0
        cache = local_empty_cache(nTargets, nSources, target_sites, source_sites, ...
            target_pos, source_pos, target_alpha, source_alpha, ...
            alphaEwald, rcut, Lmin, exclude_self, use_thole);
        return;
    end

    cutoff2 = rcut^2;

    ti_chunks = {};
    sj_chunks = {};
    dr_chunks = {};
    r2_chunks = {};

    nChunk = 0;

    for b0 = 1:target_block_size:nTargets
        b1 = min(b0 + target_block_size - 1, nTargets);
        blockIdx = (b0:b1).';
        nb = numel(blockIdx);

        Xt = target_pos(blockIdx, :);

        % Cartesian displacement from source to target:
        %   dr0 = r_target - r_source
        dx = Xt(:,1) - source_pos(:,1).';
        dy = Xt(:,2) - source_pos(:,2).';
        dz = Xt(:,3) - source_pos(:,3).';

        dr0 = [dx(:), dy(:), dz(:)];

        % Minimum-image wrap using column-vector lattice convention.
        frac = (invHcol * dr0.').';
        frac = frac - round(frac);
        drMin = (Hcol * frac.').';

        r2 = sum(drMin.^2, 2);

        keep = (r2 <= cutoff2);

        if exclude_self
            sameGlobal = (target_sites(blockIdx) == source_sites.');
            keep = keep & ~sameGlobal(:);
        end

        % Guard against exact zero distances that are not explicitly excluded.
        keep = keep & (r2 > 0);

        if ~any(keep)
            continue;
        end

        [localTargetSub, localSourceSub] = ind2sub([nb, nSources], find(keep));

        nChunk = nChunk + 1;
        ti_chunks{nChunk,1} = blockIdx(localTargetSub);
        sj_chunks{nChunk,1} = localSourceSub;
        dr_chunks{nChunk,1} = drMin(keep, :);
        r2_chunks{nChunk,1} = r2(keep);
    end

    if nChunk == 0
        ti = zeros(0,1);
        sj = zeros(0,1);
        dr = zeros(0,3);
        r2_bare = zeros(0,1);
    else
        ti = vertcat(ti_chunks{:});
        sj = vertcat(sj_chunks{:});
        dr = vertcat(dr_chunks{:});
        r2_bare = vertcat(r2_chunks{:});
    end

    r_bare = sqrt(r2_bare);

    if any(r_bare <= 0)
        error('geom:build_target_source_field_cache_periodic:NonPositiveDistance', ...
            'Encountered non-positive target-source distance.');
    end

    % ---------------------------------------------------------------------
    % Periodic Ewald real-space charge-field coefficient.
    %
    % For r_vec = r_target - r_source_image:
    %
    %   E_real = q * B_alpha(r) * r_vec
    %
    % where:
    %
    %   B_alpha(r) = erfc(alpha r)/r^3
    %              + (2 alpha/sqrt(pi)) exp(-alpha^2 r^2)/r^2.
    % ---------------------------------------------------------------------
    if isempty(r_bare)
        B_ewald = zeros(0,1);
        thole_delta = zeros(0,1);
        inv_r3_eff = zeros(0,1);
    else
        invR2 = 1 ./ r2_bare;
        invR  = 1 ./ r_bare;
        invR3 = invR .* invR2;

        alpha2 = alphaEwald^2;
        twoAlphaOverSqrtPi = 2 * alphaEwald / sqrt(pi);

        erfcar = erfc(alphaEwald * r_bare);
        exp_ar2 = exp(-alpha2 * r2_bare);

        B_ewald = erfcar .* invR3 + twoAlphaOverSqrtPi .* exp_ar2 .* invR2;

        thole_delta = zeros(size(B_ewald));

        if use_thole
            alpha_i = target_alpha(ti);
            alpha_j = source_alpha(sj);

            thole_delta = local_thole_charge_correction( ...
                r_bare, alpha_i, alpha_j, thole_a);

            inv_r3_eff = B_ewald + thole_delta;
        else
            inv_r3_eff = B_ewald;
        end
    end

    % Sort by target row to form CSR.
    [ti, order] = sort(ti);
    sj = sj(order);
    dr = dr(order, :);
    r_bare = r_bare(order);
    r2_bare = r2_bare(order);
    inv_r3_eff = inv_r3_eff(order);
    B_ewald = B_ewald(order);
    thole_delta = thole_delta(order);

    nEntries = numel(ti);
    counts = accumarray(ti, 1, [nTargets, 1], @sum, 0);
    row_ptr = zeros(nTargets + 1, 1);
    row_ptr(1) = 1;
    row_ptr(2:end) = 1 + cumsum(counts);

    cache = struct();
    cache.mode = 'periodic_target_source_field';
    cache.nTargets = nTargets;
    cache.nSources = nSources;
    cache.target_sites = target_sites;
    cache.source_sites = source_sites;
    cache.target_pos = target_pos;
    cache.source_pos = source_pos;
    cache.target_alpha = target_alpha;
    cache.source_alpha = source_alpha;
    cache.row_ptr = row_ptr;
    cache.col_idx = sj;
    cache.dr = dr;
    cache.r_bare = r_bare;
    cache.r2_bare = r2_bare;
    cache.inv_r3_eff = inv_r3_eff;
    cache.B_ewald = B_ewald;
    cache.thole_delta = thole_delta;
    cache.nEntries = nEntries;
    cache.alpha = alphaEwald;
    cache.rcut = rcut;
    cache.Lmin = Lmin;
    cache.exclude_self = exclude_self;
    cache.use_thole = use_thole;
    cache.lattice_rows = Hrow;
    cache.lattice_columns = Hcol;

    if profile
        fprintf('  nEntries = %d\n', nEntries);
        fprintf('  ||B_ewald||_2     = %.16e\n', norm(B_ewald));
        fprintf('  ||thole_delta||_2 = %.16e\n', norm(thole_delta));
        fprintf('  total time = %.6f s\n', toc(tStart));
    end
end

function cache = local_empty_cache(nTargets, nSources, target_sites, source_sites, ...
    target_pos, source_pos, target_alpha, source_alpha, ...
    alphaEwald, rcut, Lmin, exclude_self, use_thole)

    cache = struct();
    cache.mode = 'periodic_target_source_field';
    cache.nTargets = nTargets;
    cache.nSources = nSources;
    cache.target_sites = target_sites;
    cache.source_sites = source_sites;
    cache.target_pos = target_pos;
    cache.source_pos = source_pos;
    cache.target_alpha = target_alpha;
    cache.source_alpha = source_alpha;
    cache.row_ptr = ones(nTargets + 1, 1);
    cache.col_idx = zeros(0, 1);
    cache.dr = zeros(0, 3);
    cache.r_bare = zeros(0, 1);
    cache.r2_bare = zeros(0, 1);
    cache.inv_r3_eff = zeros(0, 1);
    cache.B_ewald = zeros(0, 1);
    cache.thole_delta = zeros(0, 1);
    cache.nEntries = 0;
    cache.alpha = alphaEwald;
    cache.rcut = rcut;
    cache.Lmin = Lmin;
    cache.exclude_self = exclude_self;
    cache.use_thole = use_thole;
end

function delta_inv_r3 = local_thole_charge_correction(r, alpha_i, alpha_j, thole_a)
    % Short-range correction to the Ewald Coulomb real-space charge field:
    %
    %   E = q * [B_ewald(r) + (f_thole(r)-1)/r^3] * r_vec
    %
    % where
    %   f_thole(r) = 1 - exp(-a u^3)
    %   u = r / (alpha_i alpha_j)^(1/6)
    %
    % Therefore this helper returns:
    %
    %   delta_inv_r3 = (f_thole - 1) / r^3
    %
    % If either polarizability is zero, no Thole correction is applied.

    r = r(:);
    alpha_i = alpha_i(:);
    alpha_j = alpha_j(:);

    if numel(alpha_i) ~= numel(r) || numel(alpha_j) ~= numel(r)
        error('geom:build_target_source_field_cache_periodic:BadTholeAlphaSize', ...
            'alpha_i and alpha_j must match r.');
    end

    if any(alpha_i < 0) || any(alpha_j < 0)
        error('geom:build_target_source_field_cache_periodic:NegativeAlpha', ...
            'alpha_i and alpha_j must be nonnegative.');
    end

    if isscalar(thole_a)
        if thole_a < 0
            error('geom:build_target_source_field_cache_periodic:NegativeTholeA', ...
                'thole_a must be nonnegative.');
        end
        a = thole_a * ones(size(r));
    else
        a = thole_a(:);
        if numel(a) ~= numel(r)
            error('geom:build_target_source_field_cache_periodic:VectorTholeASize', ...
                'Non-scalar thole_a must match the pair-array length.');
        end
        if any(a < 0)
            error('geom:build_target_source_field_cache_periodic:NegativeTholeA', ...
                'thole_a entries must be nonnegative.');
        end
    end

    delta_inv_r3 = zeros(size(r));

    mask = (a ~= 0) & (alpha_i ~= 0) & (alpha_j ~= 0);
    if ~any(mask)
        return;
    end

    invR3 = 1 ./ (r(mask).^3);

    u = r(mask) ./ (alpha_i(mask) .* alpha_j(mask)).^(1/6);
    f = 1 - exp(-a(mask) .* u.^3);

    delta_inv_r3(mask) = (f - 1) .* invR3;
end

function x = local_require_field(s, name)
    if ~isfield(s, name) || isempty(s.(name))
        error('geom:build_target_source_field_cache_periodic:MissingField', ...
            'sys.%s is required.', name);
    end
    x = s.(name);
end

function x = local_require_struct_field(s, name)
    if ~isfield(s, name) || isempty(s.(name))
        error('geom:build_target_source_field_cache_periodic:MissingEwaldField', ...
            'ewaldParams.%s is required.', name);
    end
    x = s.(name);
end

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('geom:build_target_source_field_cache_periodic:MissingLattice', ...
            'sys.super_lattice or sys.lattice is required.');
    end

    validateattributes(H, {'double'}, {'size',[3,3],'real','finite'}, ...
        mfilename, 'lattice');
end