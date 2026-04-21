function cache = build_periodic_realspace_cache(sys, problem, ewaldParams, scfParams)
%BUILD_PERIODIC_REALSPACE_CACHE Build periodic real-space Ewald interaction cache.
%
% Chunked one-pass version:
%   - avoids separate count pass
%   - avoids AGROW inside the hot inner loop
%   - accumulates one chunk per outer aa and concatenates once
%
% Output format matches the previous implementation.

    narginchk(3, 4);

    if nargin < 4 || isempty(scfParams)
        scfParams = struct();
    end

    validateattributes(problem, {'struct'}, {'scalar'}, mfilename, 'problem', 2);
    validateattributes(ewaldParams, {'struct'}, {'scalar'}, mfilename, 'ewaldParams', 3);

    verbose = false;
    if isfield(scfParams, 'verbose') && ~isempty(scfParams.verbose)
        verbose = logical(scfParams.verbose);
    end

    tTotal = tic;

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('geom:build_periodic_realspace_cache:MissingSitePos', ...
            'sys.site_pos is required and may not be empty.');
    end
    pos = sys.site_pos;
    validateattributes(pos, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
        mfilename, 'sys.site_pos');

    nSites = size(pos, 1);

    if ~isfield(problem, 'polMask') || isempty(problem.polMask)
        error('geom:build_periodic_realspace_cache:MissingPolMask', ...
            'problem.polMask is required.');
    end
    siteMask = logical(problem.polMask(:));
    if numel(siteMask) ~= nSites
        error('geom:build_periodic_realspace_cache:BadPolMask', ...
            'problem.polMask must have length equal to size(sys.site_pos,1).');
    end
    activeSites = find(siteMask);
    nActive = numel(activeSites);

    if ~isfield(ewaldParams, 'alpha') || isempty(ewaldParams.alpha)
        error('geom:build_periodic_realspace_cache:MissingAlpha', ...
            'ewaldParams.alpha is required.');
    end
    alpha = ewaldParams.alpha;
    validateattributes(alpha, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'ewaldParams.alpha');

    if ~isfield(ewaldParams, 'rcut') || isempty(ewaldParams.rcut)
        error('geom:build_periodic_realspace_cache:MissingRcut', ...
            'ewaldParams.rcut is required.');
    end
    rcut = ewaldParams.rcut;
    validateattributes(rcut, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'ewaldParams.rcut');

    use_thole = true;
    if isfield(scfParams, 'use_thole') && ~isempty(scfParams.use_thole)
        use_thole = logical(scfParams.use_thole);
    end

    if use_thole
        if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
            error('geom:build_periodic_realspace_cache:MissingSiteAlpha', ...
                'sys.site_alpha is required when use_thole = true.');
        end
        if ~isfield(sys, 'thole_a') || isempty(sys.thole_a)
            error('geom:build_periodic_realspace_cache:MissingTholeA', ...
                'sys.thole_a is required when use_thole = true.');
        end
        alpha_site = sys.site_alpha(:);
        if numel(alpha_site) ~= nSites
            error('geom:build_periodic_realspace_cache:BadSiteAlpha', ...
                'sys.site_alpha must have length equal to size(sys.site_pos,1).');
        end
        thole_a = sys.thole_a;
    else
        alpha_site = [];
        thole_a = [];
    end

    H = local_get_direct_lattice(sys);
    a = H(:,1);
    b = H(:,2);
    c = H(:,3);

    nxmax = ceil(rcut / norm(a)) + 1;
    nymax = ceil(rcut / norm(b)) + 1;
    nzmax = ceil(rcut / norm(c)) + 1;

    % Enumerate all image shifts once.
    tEnum = tic;

    nAllocR = (2*nxmax + 1) * (2*nymax + 1) * (2*nzmax + 1);
    imageN = zeros(nAllocR, 3);
    Rshifts = zeros(nAllocR, 3);
    isZeroShift = false(nAllocR, 1);

    nR = 0;
    for nx = -nxmax:nxmax
        for ny = -nymax:nymax
            for nz = -nzmax:nzmax
                nR = nR + 1;
                nvec = [nx; ny; nz];
                imageN(nR, :) = nvec.';
                Rshifts(nR, :) = (H * nvec).';
                isZeroShift(nR) = (nx == 0 && ny == 0 && nz == 0);
            end
        end
    end

    imageN = imageN(1:nR, :);
    Rshifts = Rshifts(1:nR, :);
    isZeroShift = isZeroShift(1:nR);

    timeEnum = toc(tEnum);

    rcut2 = rcut^2;
    alpha2 = alpha^2;
    twoAlphaOverSqrtPi = 2 * alpha / sqrt(pi);

    if verbose
        fprintf(['build_periodic_realspace_cache: nActive=%d | rcut=%.6f | ' ...
                 'image bounds=[%d %d %d] | nImageShifts=%d\n'], ...
            nActive, rcut, nxmax, nymax, nzmax, nR);
        fprintf('  image enumeration done in %.3f s\n', timeEnum);
    end

    % ------------------------------------------------------------------
    % Distinct-site accumulation: one chunk per outer aa

    tDistinct = tic;

    nChunk = max(nActive - 1, 1);
    pair_i_cells      = cell(nChunk, 1);
    pair_j_cells      = cell(nChunk, 1);
    image_n_cells     = cell(nChunk, 1);
    dr_cells          = cell(nChunk, 1);
    r_bare_cells      = cell(nChunk, 1);
    r2_bare_cells     = cell(nChunk, 1);
    coeff_iso_cells   = cell(nChunk, 1);
    coeff_dyad_cells  = cell(nChunk, 1);

    nDistinctKept = 0;

    % Heuristic initial capacity per aa chunk.
    % Large enough to reduce reallocations, small enough not to explode memory.
    initCap = max(4096, 8 * nR);

    for aa = 1:(nActive - 1)
        i = activeSites(aa);
        ri = pos(i, :);

        cap = initCap;
        cursor = 0;

        pair_i_blk     = zeros(cap, 1);
        pair_j_blk     = zeros(cap, 1);
        image_n_blk    = zeros(cap, 3);
        dr_blk         = zeros(cap, 3);
        r_bare_blk     = zeros(cap, 1);
        r2_bare_blk    = zeros(cap, 1);
        coeff_iso_blk  = zeros(cap, 1);
        coeff_dyad_blk = zeros(cap, 1);

        for bb = (aa + 1):nActive
            j = activeSites(bb);
            rij0 = pos(j, :) - ri;

            xvec = rij0 + Rshifts;
            x2 = sum(xvec.^2, 2);
            keep = (x2 > 0) & (x2 <= rcut2);

            if ~any(keep)
                continue;
            end

            x2keep = x2(keep);
            xkeep = sqrt(x2keep);
            xveckeep = xvec(keep, :);
            imageKeep = imageN(keep, :);

            erfcax = erfc(alpha * xkeep);
            expax2 = exp(-alpha2 * x2keep);

            invx2 = 1 ./ x2keep;
            invx = 1 ./ xkeep;
            invx3 = invx .* invx2;
            invx5 = invx3 .* invx2;
            invx4 = invx2.^2;

            B = erfcax .* invx3 + twoAlphaOverSqrtPi * expax2 .* invx2;
            C = 3 .* erfcax .* invx5 + twoAlphaOverSqrtPi .* ...
                (2 .* alpha2 .* invx2 + 3 .* invx4) .* expax2;

            coeff_iso = -B;
            coeff_dyad = +C;

            if use_thole
                tf = local_thole_factors_vectorized(xkeep, alpha_site(i), alpha_site(j), thole_a);
                coeff_iso = coeff_iso - tf.l3 .* invx3;
                coeff_dyad = coeff_dyad + 3 .* tf.l5 .* invx5;
            end

            m = numel(xkeep);
            need = cursor + m;

            if need > cap
                newCap = cap;
                while need > newCap
                    newCap = 2 * newCap;
                end

                pair_i_blk(newCap, 1) = 0;
                pair_j_blk(newCap, 1) = 0;
                image_n_blk(newCap, 3) = 0;
                dr_blk(newCap, 3) = 0;
                r_bare_blk(newCap, 1) = 0;
                r2_bare_blk(newCap, 1) = 0;
                coeff_iso_blk(newCap, 1) = 0;
                coeff_dyad_blk(newCap, 1) = 0;

                cap = newCap;
            end

            idx = (cursor + 1):(cursor + m);

            pair_i_blk(idx) = i;
            pair_j_blk(idx) = j;
            image_n_blk(idx, :) = imageKeep;
            dr_blk(idx, :) = xveckeep;
            r_bare_blk(idx) = xkeep;
            r2_bare_blk(idx) = x2keep;
            coeff_iso_blk(idx) = coeff_iso;
            coeff_dyad_blk(idx) = coeff_dyad;

            cursor = need;
            nDistinctKept = nDistinctKept + m;
        end

        pair_i_cells{aa}     = pair_i_blk(1:cursor);
        pair_j_cells{aa}     = pair_j_blk(1:cursor);
        image_n_cells{aa}    = image_n_blk(1:cursor, :);
        dr_cells{aa}         = dr_blk(1:cursor, :);
        r_bare_cells{aa}     = r_bare_blk(1:cursor);
        r2_bare_cells{aa}    = r2_bare_blk(1:cursor);
        coeff_iso_cells{aa}  = coeff_iso_blk(1:cursor);
        coeff_dyad_cells{aa} = coeff_dyad_blk(1:cursor);
    end

    timeDistinct = toc(tDistinct);

    % ------------------------------------------------------------------
    % Self-image accumulation

    tSelf = tic;

    keepSelf = ~isZeroShift;
    Rself = Rshifts(keepSelf, :);
    imageSelf = imageN(keepSelf, :);

    pair_i_self = zeros(0,1);
    pair_j_self = zeros(0,1);
    image_n_self = zeros(0,3);
    dr_self = zeros(0,3);
    r_bare_self = zeros(0,1);
    r2_bare_self = zeros(0,1);
    coeff_iso_self_all = zeros(0,1);
    coeff_dyad_self_all = zeros(0,1);

    nSelfKept = 0;
    nSelfPerSite = 0;

    if ~isempty(Rself)
        x2 = sum(Rself.^2, 2);
        keep = (x2 > 0) & (x2 <= rcut2);

        if any(keep)
            Rself = Rself(keep, :);
            imageSelf = imageSelf(keep, :);
            x2 = x2(keep);
            x = sqrt(x2);
            nSelfPerSite = numel(x);

            erfcax = erfc(alpha * x);
            expax2 = exp(-alpha2 * x2);

            invx2 = 1 ./ x2;
            invx = 1 ./ x;
            invx3 = invx .* invx2;
            invx5 = invx3 .* invx2;
            invx4 = invx2.^2;

            B = erfcax .* invx3 + twoAlphaOverSqrtPi * expax2 .* invx2;
            C = 3 .* erfcax .* invx5 + twoAlphaOverSqrtPi .* ...
                (2 .* alpha2 .* invx2 + 3 .* invx4) .* expax2;

            coeff_iso_base = -B;
            coeff_dyad_base = +C;

            nSelfTotal = nActive * nSelfPerSite;
            pair_i_self = zeros(nSelfTotal, 1);
            pair_j_self = zeros(nSelfTotal, 1);
            image_n_self = zeros(nSelfTotal, 3);
            dr_self = zeros(nSelfTotal, 3);
            r_bare_self = zeros(nSelfTotal, 1);
            r2_bare_self = zeros(nSelfTotal, 1);
            coeff_iso_self_all = zeros(nSelfTotal, 1);
            coeff_dyad_self_all = zeros(nSelfTotal, 1);

            cursor = 0;
            for aa = 1:nActive
                i = activeSites(aa);

                coeff_iso = coeff_iso_base;
                coeff_dyad = coeff_dyad_base;

                if use_thole
                    tf = local_thole_factors_vectorized(x, alpha_site(i), alpha_site(i), thole_a);
                    coeff_iso = coeff_iso - tf.l3 .* invx3;
                    coeff_dyad = coeff_dyad + 3 .* tf.l5 .* invx5;
                end

                m = nSelfPerSite;
                idx = (cursor + 1):(cursor + m);

                pair_i_self(idx) = i;
                pair_j_self(idx) = i;
                image_n_self(idx, :) = imageSelf;
                dr_self(idx, :) = Rself;
                r_bare_self(idx) = x;
                r2_bare_self(idx) = x2;
                coeff_iso_self_all(idx) = coeff_iso;
                coeff_dyad_self_all(idx) = coeff_dyad;

                cursor = cursor + m;
            end

            nSelfKept = nSelfTotal;
        end
    end

    timeSelf = toc(tSelf);

    % ------------------------------------------------------------------
    % Concatenate once

    tConcat = tic;

    if nActive > 1
        pair_i = vertcat(pair_i_cells{:});
        pair_j = vertcat(pair_j_cells{:});
        image_n = vertcat(image_n_cells{:});
        dr_all = vertcat(dr_cells{:});
        r_bare_all = vertcat(r_bare_cells{:});
        r2_bare_all = vertcat(r2_bare_cells{:});
        coeff_iso_all = vertcat(coeff_iso_cells{:});
        coeff_dyad_all = vertcat(coeff_dyad_cells{:});
    else
        pair_i = zeros(0,1);
        pair_j = zeros(0,1);
        image_n = zeros(0,3);
        dr_all = zeros(0,3);
        r_bare_all = zeros(0,1);
        r2_bare_all = zeros(0,1);
        coeff_iso_all = zeros(0,1);
        coeff_dyad_all = zeros(0,1);
    end

    if nSelfKept > 0
        pair_i = [pair_i; pair_i_self];
        pair_j = [pair_j; pair_j_self];
        image_n = [image_n; image_n_self];
        dr_all = [dr_all; dr_self];
        r_bare_all = [r_bare_all; r_bare_self];
        r2_bare_all = [r2_bare_all; r2_bare_self];
        coeff_iso_all = [coeff_iso_all; coeff_iso_self_all];
        coeff_dyad_all = [coeff_dyad_all; coeff_dyad_self_all];
    end

    timeConcat = toc(tConcat);

    nKeep = numel(pair_i);

    cache = struct();
    cache.mode = 'periodic_realspace';
    cache.nSites = nSites;
    cache.site_mask = siteMask;
    cache.site_idx = activeSites;
    cache.activeSites = activeSites;

    cache.pair_i = pair_i;
    cache.pair_j = pair_j;
    cache.image_n = image_n;
    cache.dr = dr_all;

    cache.r_bare = r_bare_all;
    cache.r2_bare = r2_bare_all;

    cache.coeff_iso = coeff_iso_all;
    cache.coeff_dyad = coeff_dyad_all;

    cache.alpha = alpha;
    cache.rcut = rcut;
    cache.use_thole = use_thole;

    cache.real_image_bounds = [nxmax, nymax, nzmax];
    cache.nImageShifts = nR;
    cache.nInteractions = nKeep;
    cache.includes_self_image = true;

    cache.timing = struct();
    cache.timing.image_enum_s     = timeEnum;
    cache.timing.distinct_fill_s  = timeDistinct;
    cache.timing.self_fill_s      = timeSelf;
    cache.timing.concat_s         = timeConcat;
    cache.timing.total_s          = toc(tTotal);
    cache.timing.nDistinctKept    = nDistinctKept;
    cache.timing.nSelfKept        = nSelfKept;
    cache.timing.nSelfPerSite     = nSelfPerSite;

    if verbose
        fprintf('  distinct fill done in %.3f s | kept=%d\n', timeDistinct, nDistinctKept);
        fprintf('  self-image fill done in %.3f s | nSelfPerSite=%d | kept=%d\n', ...
            timeSelf, nSelfPerSite, nSelfKept);
        fprintf('  concatenate done in %.3f s\n', timeConcat);
        fprintf('build_periodic_realspace_cache: total time = %.3f s\n', cache.timing.total_s);
    end
end

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('geom:build_periodic_realspace_cache:MissingLattice', ...
            'sys.super_lattice or sys.lattice is required for periodic real-space cache.');
    end

    validateattributes(H, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
        mfilename, 'lattice');
end

function tf = local_thole_factors_vectorized(r, alpha_i, alpha_j, thole_a)
    if alpha_i < 0 || alpha_j < 0
        error('geom:build_periodic_realspace_cache:NegativeAlpha', ...
            'alpha_i and alpha_j must be nonnegative.');
    end
    if thole_a < 0
        error('geom:build_periodic_realspace_cache:NegativeTholeA', ...
            'thole_a must be nonnegative.');
    end

    r = r(:);

    tf = struct();
    tf.f3 = ones(size(r));
    tf.f5 = ones(size(r));
    tf.l3 = zeros(size(r));
    tf.l5 = zeros(size(r));

    if thole_a == 0
        return;
    end

    if alpha_i == 0 || alpha_j == 0
        return;
    end

    u = r ./ (alpha_i * alpha_j)^(1/6);
    au3 = thole_a * u.^3;

    tf.f3 = 1 - exp(-au3);
    tf.f5 = 1 - (1 + au3) .* exp(-au3);
    tf.l3 = tf.f3 - 1;
    tf.l5 = tf.f5 - 1;
end