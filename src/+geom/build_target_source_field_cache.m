function cache = build_target_source_field_cache(sys, opts)
%BUILD_TARGET_SOURCE_FIELD_CACHE Build asymmetric target-source cache for
%nonperiodic external field from charges.
%
% cache = geom.build_target_source_field_cache(sys, opts)
%
% Inputs
%   sys   canonical polarization system in atomic units
%
%   opts  struct with fields:
%     .target_mask   logical nSites x 1, required
%     .source_mask   logical nSites x 1, required
%     .exclude_self  logical scalar, default true
%     .rcut          scalar cutoff in bohr, default inf
%     .use_thole     logical scalar, default false
%     .softening     scalar, default 0.0
%
% Output
%   cache struct with fields:
%     .kind              = 'target_source_field_cache'
%     .mode              = 'nonperiodic'
%     .nSites
%     .nPairs
%     .pair_target_full
%     .pair_source_full
%     .dr                = r_source - r_target
%     .r2_bare
%     .inv_r3_bare
%     .thole_f3          optional
%     .rcut
%     .exclude_self

    narginchk(2, 2);
    io.assert_atomic_units(sys);

    if ~isfield(opts, 'target_mask') || isempty(opts.target_mask)
        error('geom:build_target_source_field_cache:MissingTargetMask', ...
            'opts.target_mask is required.');
    end
    if ~isfield(opts, 'source_mask') || isempty(opts.source_mask)
        error('geom:build_target_source_field_cache:MissingSourceMask', ...
            'opts.source_mask is required.');
    end

    nSites = sys.n_sites;
    target_mask = logical(opts.target_mask(:));
    source_mask = logical(opts.source_mask(:));

    if numel(target_mask) ~= nSites || numel(source_mask) ~= nSites
        error('geom:build_target_source_field_cache:MaskSizeMismatch', ...
            'target_mask and source_mask must have length sys.n_sites.');
    end

    exclude_self = true;
    if isfield(opts, 'exclude_self') && ~isempty(opts.exclude_self)
        exclude_self = logical(opts.exclude_self);
    end

    rcut = inf;
    if isfield(opts, 'rcut') && ~isempty(opts.rcut)
        rcut = opts.rcut;
    end
    rcut2 = rcut^2;

    useThole = false;
    if isfield(opts, 'use_thole') && ~isempty(opts.use_thole)
        useThole = logical(opts.use_thole);
    end

    softening = 0.0;
    if isfield(opts, 'softening') && ~isempty(opts.softening)
        softening = opts.softening;
    end

    pos = sys.site_pos;

    targetIdx = find(target_mask);

    % Only charged sources can contribute to Eext, so filter them here.
    if ~isfield(sys, 'site_charge') || isempty(sys.site_charge)
        error('geom:build_target_source_field_cache:MissingSiteCharge', ...
            'sys.site_charge is required.');
    end
    q = sys.site_charge(:);
    sourceIdx = find(source_mask & (q ~= 0));

    nTarget = numel(targetIdx);
    nSource = numel(sourceIdx);

    if nTarget == 0 || nSource == 0
        cache = struct();
        cache.kind = 'target_source_field_cache';
        cache.mode = 'nonperiodic';
        cache.nSites = nSites;
        cache.nPairs = 0;
        cache.pair_target_full = zeros(0, 1);
        cache.pair_source_full = zeros(0, 1);
        cache.dr = zeros(0, 3);
        cache.r2_bare = zeros(0, 1);
        cache.inv_r3_bare = zeros(0, 1);
        if useThole
            cache.thole_f3 = zeros(0, 1);
        end
        cache.rcut = rcut;
        cache.exclude_self = exclude_self;
        return;
    end

    targetPos = pos(targetIdx, :);

    pairTargetBlocks = cell(nSource, 1);
    pairSourceBlocks = cell(nSource, 1);
    drBlocks         = cell(nSource, 1);
    r2Blocks         = cell(nSource, 1);
    invR3Blocks      = cell(nSource, 1);
    if useThole
        f3Blocks = cell(nSource, 1);
        alphaSites = sys.site_alpha(:);
        tholeA = sys.thole_a;
    end

    nPairs = 0;

    for s = 1:nSource
        j = sourceIdx(s);
        rj = pos(j, :);

        % dr = r_source - r_target for all targets
        dr = rj - targetPos;                    % nTarget x 3
        r2_bare = sum(dr.^2, 2);               % nTarget x 1

        keep = true(nTarget, 1);

        if exclude_self
            keep = keep & (targetIdx ~= j);
        end

        if isfinite(rcut)
            keep = keep & (r2_bare <= rcut2);
        end

        if ~any(keep)
            pairTargetBlocks{s} = zeros(0, 1);
            pairSourceBlocks{s} = zeros(0, 1);
            drBlocks{s}         = zeros(0, 3);
            r2Blocks{s}         = zeros(0, 1);
            invR3Blocks{s}      = zeros(0, 1);
            if useThole
                f3Blocks{s} = zeros(0, 1);
            end
            continue;
        end

        dr_keep = dr(keep, :);
        r2_keep = r2_bare(keep);
        tgt_keep = targetIdx(keep);

        nKeep = numel(tgt_keep);
        nPairs = nPairs + nKeep;

        pairTargetBlocks{s} = tgt_keep;
        pairSourceBlocks{s} = repmat(j, nKeep, 1);
        drBlocks{s}         = dr_keep;
        r2Blocks{s}         = r2_keep;

        if softening == 0
            r_keep = sqrt(r2_keep);
            invR3Blocks{s} = 1 ./ (r_keep.^3);
        else
            invR3Blocks{s} = NaN(nKeep, 1);
            r_keep = sqrt(r2_keep);
        end

        if useThole
            % Scalar Thole helper applied over the kept target list.
            alpha_j = alphaSites(j);
            f3 = zeros(nKeep, 1);
            for k = 1:nKeep
                i = tgt_keep(k);
                tf = thole.thole_f3f5_factors(r_keep(k), alphaSites(i), alpha_j, tholeA);
                f3(k) = tf.f3;
            end
            f3Blocks{s} = f3;
        end
    end

    cache = struct();
    cache.kind = 'target_source_field_cache';
    cache.mode = 'nonperiodic';
    cache.nSites = nSites;
    cache.nPairs = nPairs;

    if nPairs == 0
        cache.pair_target_full = zeros(0, 1);
        cache.pair_source_full = zeros(0, 1);
        cache.dr = zeros(0, 3);
        cache.r2_bare = zeros(0, 1);
        cache.inv_r3_bare = zeros(0, 1);
        if useThole
            cache.thole_f3 = zeros(0, 1);
        end
    else
        cache.pair_target_full = vertcat(pairTargetBlocks{:});
        cache.pair_source_full = vertcat(pairSourceBlocks{:});
        cache.dr               = vertcat(drBlocks{:});
        cache.r2_bare          = vertcat(r2Blocks{:});
        cache.inv_r3_bare      = vertcat(invR3Blocks{:});
        if useThole
            cache.thole_f3 = vertcat(f3Blocks{:});
        end
    end

    cache.rcut = rcut;
    cache.exclude_self = exclude_self;
end