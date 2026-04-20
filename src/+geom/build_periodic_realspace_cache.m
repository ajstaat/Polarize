function cache = build_periodic_realspace_cache(sys, problem, ewaldParams, scfParams)
%BUILD_PERIODIC_REALSPACE_CACHE Build periodic real-space Ewald interaction cache.
%
% cache = geom.build_periodic_realspace_cache(sys, problem, ewaldParams)
% cache = geom.build_periodic_realspace_cache(sys, problem, ewaldParams, scfParams)
%
% Inputs
%   sys         canonical polarization-system struct in atomic units
%   problem     struct from thole.prepare_scf_problem(...)
%   ewaldParams struct with fields:
%       .alpha     Ewald screening parameter
%       .rcut      real-space cutoff
%
%   scfParams   optional struct with fields:
%       .use_thole   logical, default true
%
% Output
%   cache       struct with periodic real-space pair-image interactions
%
% Fields
%   .mode               'periodic_realspace'
%   .nSites
%   .site_mask
%   .site_idx
%   .activeSites
%   .pair_i             full-site target index
%   .pair_j             full-site source index
%   .image_n            integer image coefficients [nx ny nz]
%   .dr                 xvec = r_j + R - r_i
%   .r_bare
%   .r2_bare
%   .coeff_iso          scalar coefficient multiplying mu_j
%   .coeff_dyad         scalar coefficient multiplying xvec*(xvec·mu_j)
%   .alpha
%   .rcut
%   .use_thole
%   .real_image_bounds
%   .nImageShifts
%   .nInteractions
%   .includes_self_image
%
% Notes
%   - This is the periodic real-space analog of the nonperiodic pair cache.
%   - The tensor action is represented as
%         E_{i<-j} = coeff_iso * mu_j + coeff_dyad * xvec * (xvec · mu_j)
%     where xvec = r_j + R - r_i.
%   - coeff_iso / coeff_dyad already include the optional real-space Thole
%     correction when scfParams.use_thole is true.
%   - Distinct-site interactions are stored once per unordered pair-image
%     object (i < j), plus self-image entries with i == j and R ~= 0.
%   - Central self (i == j, R == 0) is excluded.
%   - This implementation uses a two-pass strategy:
%       pass 1: count surviving interactions
%       pass 2: allocate exactly and fill
%     to avoid pathological worst-case memory allocation.

    narginchk(3, 4);

    if nargin < 4 || isempty(scfParams)
        scfParams = struct();
    end

    validateattributes(problem, {'struct'}, {'scalar'}, mfilename, 'problem', 2);
    validateattributes(ewaldParams, {'struct'}, {'scalar'}, mfilename, 'ewaldParams', 3);

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

    rcut2 = rcut^2;

    alpha2 = alpha^2;
    twoAlphaOverSqrtPi = 2 * alpha / sqrt(pi);

    % ------------------------------------------------------------------
    % PASS 1: count exact number of surviving interactions

    nKeep = 0;

    % Distinct-site interactions: i < j
    for aa = 1:(nActive - 1)
        i = activeSites(aa);
        ri = pos(i, :);

        for bb = (aa + 1):nActive
            j = activeSites(bb);
            rij0 = pos(j, :) - ri;      % r_j - r_i

            xvec = rij0 + Rshifts;      % nR x 3
            x2 = sum(xvec.^2, 2);

            keep = (x2 > 0) & (x2 <= rcut2);
            nKeep = nKeep + nnz(keep);
        end
    end

    % Self-image interactions: i == j, R ~= 0 only
    keepSelf = ~isZeroShift;
    Rself = Rshifts(keepSelf, :);
    imageSelf = imageN(keepSelf, :);

    if ~isempty(Rself)
        x2 = sum(Rself.^2, 2);
        keep = (x2 > 0) & (x2 <= rcut2);
        nSelfPerSite = nnz(keep);
        nKeep = nKeep + nActive * nSelfPerSite;
    end

    % Allocate exactly
    pair_i = zeros(nKeep, 1);
    pair_j = zeros(nKeep, 1);
    image_n = zeros(nKeep, 3);
    dr_all = zeros(nKeep, 3);
    r_bare_all = zeros(nKeep, 1);
    r2_bare_all = zeros(nKeep, 1);
    coeff_iso_all = zeros(nKeep, 1);
    coeff_dyad_all = zeros(nKeep, 1);

    % ------------------------------------------------------------------
    % PASS 2: fill arrays

    nFill = 0;

    % Distinct-site interactions: i < j, all image shifts inside cutoff.
    for aa = 1:(nActive - 1)
        i = activeSites(aa);
        ri = pos(i, :);

        for bb = (aa + 1):nActive
            j = activeSites(bb);
            rij0 = pos(j, :) - ri;      % r_j - r_i

            xvec = rij0 + Rshifts;      % nR x 3
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

            % Field-tensor convention matching thole.dipole_tensor_block:
            %   T = -B I + C xx^T
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
            idx = (nFill + 1):(nFill + m);

            pair_i(idx) = i;
            pair_j(idx) = j;
            image_n(idx, :) = imageKeep;
            dr_all(idx, :) = xveckeep;
            r_bare_all(idx) = xkeep;
            r2_bare_all(idx) = x2keep;
            coeff_iso_all(idx) = coeff_iso;
            coeff_dyad_all(idx) = coeff_dyad;

            nFill = nFill + m;
        end
    end

    % Self-image interactions: i == j, R ~= 0 only.
    keepSelf = ~isZeroShift;
    Rself = Rshifts(keepSelf, :);
    imageSelf = imageN(keepSelf, :);

    if ~isempty(Rself)
        x2 = sum(Rself.^2, 2);
        keep = (x2 > 0) & (x2 <= rcut2);

        if any(keep)
            Rself = Rself(keep, :);
            imageSelf = imageSelf(keep, :);
            x2 = x2(keep);
            x = sqrt(x2);

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

            coeff_iso_self = -B;
            coeff_dyad_self = +C;

            for aa = 1:nActive
                i = activeSites(aa);

                coeff_iso = coeff_iso_self;
                coeff_dyad = coeff_dyad_self;

                if use_thole
                    tf = local_thole_factors_vectorized(x, alpha_site(i), alpha_site(i), thole_a);
                    coeff_iso = coeff_iso - tf.l3 .* invx3;
                    coeff_dyad = coeff_dyad + 3 .* tf.l5 .* invx5;
                end

                m = numel(x);
                idx = (nFill + 1):(nFill + m);

                pair_i(idx) = i;
                pair_j(idx) = i;
                image_n(idx, :) = imageSelf;
                dr_all(idx, :) = Rself;
                r_bare_all(idx) = x;
                r2_bare_all(idx) = x2;
                coeff_iso_all(idx) = coeff_iso;
                coeff_dyad_all(idx) = coeff_dyad;

                nFill = nFill + m;
            end
        end
    end

    if nFill ~= nKeep
        error('geom:build_periodic_realspace_cache:FillCountMismatch', ...
            'Two-pass fill count mismatch: counted %d, filled %d.', nKeep, nFill);
    end

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
end

function H = local_get_direct_lattice(sys)
%LOCAL_GET_DIRECT_LATTICE Extract direct lattice matrix, columns are lattice vectors.

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
%LOCAL_THOLE_FACTORS_VECTORIZED Vectorized Thole factors for one site pair and many distances.
%
% Mirrors thole.thole_f3f5_factors but for vector input r.

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