function cache = build_periodic_realspace_geom_cache(sys, problem, ewaldParams, opts)
%BUILD_PERIODIC_REALSPACE_GEOM_CACHE Build periodic real-space pair-image geometry cache.
%
% cache = geom.build_periodic_realspace_geom_cache(sys, problem, ewaldParams)
% cache = geom.build_periodic_realspace_geom_cache(sys, problem, ewaldParams, opts)
%
% Fast path:
%   - uses mex_build_periodic_realspace_geom_cache when available
%
% Fallback:
%   - preserves the uploaded MATLAB implementation structure
%
% Inputs
%   sys         canonical polarization-system struct in atomic units
%   problem     struct from thole.prepare_scf_problem(...)
%   ewaldParams struct with fields:
%       .rcut      real-space cutoff
%
%   opts        optional struct with fields:
%       .target_mask   N x 1 logical, optional
%       .source_mask   N x 1 logical, optional
%
% Output
%   cache       struct with periodic real-space pair-image geometry
%
% Fields
%   .mode               'periodic_realspace_geom'
%   .nSites
%   .target_mask
%   .source_mask
%   .target_idx
%   .source_idx
%   .activeSites        alias of target_idx for backward compatibility
%   .pair_i             full-site target index
%   .pair_j             full-site source index
%   .image_n            integer image coefficients [nx ny nz]
%   .dr                 xvec = r_j + R - r_i
%   .r_bare
%   .r2_bare
%   .rcut
%   .real_image_bounds
%   .nImageShifts
%   .nInteractions
%   .includes_self_image
%   .is_symmetric_pair_cache
%
% Notes
%   - If target/source masks are identical, distinct-site interactions are
%     stored once per unordered pair-image object (i < j), plus self-image
%     entries with i == j and R ~= 0, matching the old cache.
%   - If target/source masks differ, interactions are stored DIRECTED:
%         pair_i = target, pair_j = source
%     for all kept target-source-image combinations.
%   - Central self (i == j, R == 0) is excluded.
%   - This is geometry only; no dipole/charge coefficients are stored.

    narginchk(3, 4);

    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    validateattributes(problem, {'struct'}, {'scalar'}, mfilename, 'problem', 2);
    validateattributes(ewaldParams, {'struct'}, {'scalar'}, mfilename, 'ewaldParams', 3);

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('geom:build_periodic_realspace_geom_cache:MissingSitePos', ...
            'sys.site_pos is required and may not be empty.');
    end
    pos = sys.site_pos;
    validateattributes(pos, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
        mfilename, 'sys.site_pos');

    nSites = size(pos, 1);

    % Default target/source masks from problem.polMask, matching old behavior.
    if isfield(opts, 'target_mask') && ~isempty(opts.target_mask)
        targetMask = logical(opts.target_mask(:));
    else
        if ~isfield(problem, 'polMask') || isempty(problem.polMask)
            error('geom:build_periodic_realspace_geom_cache:MissingPolMask', ...
                'problem.polMask is required when opts.target_mask is not provided.');
        end
        targetMask = logical(problem.polMask(:));
    end

    if isfield(opts, 'source_mask') && ~isempty(opts.source_mask)
        sourceMask = logical(opts.source_mask(:));
    else
        if ~isfield(problem, 'polMask') || isempty(problem.polMask)
            error('geom:build_periodic_realspace_geom_cache:MissingPolMask', ...
                'problem.polMask is required when opts.source_mask is not provided.');
        end
        sourceMask = logical(problem.polMask(:));
    end

    if numel(targetMask) ~= nSites
        error('geom:build_periodic_realspace_geom_cache:BadTargetMask', ...
            'target_mask must have length equal to size(sys.site_pos,1).');
    end
    if numel(sourceMask) ~= nSites
        error('geom:build_periodic_realspace_geom_cache:BadSourceMask', ...
            'source_mask must have length equal to size(sys.site_pos,1).');
    end

    targetSites = find(targetMask);
    sourceSites = find(sourceMask);

    if ~isfield(ewaldParams, 'rcut') || isempty(ewaldParams.rcut)
        error('geom:build_periodic_realspace_geom_cache:MissingRcut', ...
            'ewaldParams.rcut is required.');
    end
    rcut = ewaldParams.rcut;
    validateattributes(rcut, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'ewaldParams.rcut');

    H = local_get_direct_lattice(sys);
    a = H(:,1);
    b = H(:,2);
    c = H(:,3);

    nxmax = ceil(rcut / norm(a)) + 1;
    nymax = ceil(rcut / norm(b)) + 1;
    nzmax = ceil(rcut / norm(c)) + 1;

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
    isSymmetric = isequal(targetMask, sourceMask);

    mexAvailable = (exist(['mex_build_periodic_realspace_geom_cache.' mexext], 'file') == 3) || ...
                   (exist('mex_build_periodic_realspace_geom_cache', 'file') == 3);

    if mexAvailable
        [pair_i, pair_j, image_n, dr_all, r_bare_all, r2_bare_all, nKeep] = ...
            mex_build_periodic_realspace_geom_cache( ...
                double(pos), ...
                double(targetMask), ...
                double(sourceMask), ...
                double(H), ...
                double(rcut));

        cache = struct();
        cache.mode = 'periodic_realspace_geom';
        cache.nSites = nSites;

        cache.target_mask = targetMask;
        cache.source_mask = sourceMask;
        cache.target_idx = targetSites;
        cache.source_idx = sourceSites;

        % Backward compatibility
        cache.site_mask = targetMask;
        cache.site_idx = targetSites;
        cache.activeSites = targetSites;

        cache.pair_i = pair_i;
        cache.pair_j = pair_j;
        cache.image_n = image_n;
        cache.dr = dr_all;

        cache.r_bare = r_bare_all;
        cache.r2_bare = r2_bare_all;

        cache.rcut = rcut;
        cache.real_image_bounds = [nxmax, nymax, nzmax];
        cache.nImageShifts = nR;
        cache.nInteractions = nKeep;
        cache.includes_self_image = true;
        cache.is_symmetric_pair_cache = isSymmetric;
        return;
    end

    % ------------------------------------------------------------------
    % MATLAB fallback: preserve uploaded implementation structure

    % PASS 1: count exact number of interactions
    nKeep = 0;

    if isSymmetric
        nTarget = numel(targetSites);

        for aa = 1:(nTarget - 1)
            i = targetSites(aa);
            ri = pos(i, :);

            for bb = (aa + 1):nTarget
                j = targetSites(bb);
                rij0 = pos(j, :) - ri;

                xvec = rij0 + Rshifts;
                x2 = sum(xvec.^2, 2);

                keep = (x2 > 0) & (x2 <= rcut2);
                nKeep = nKeep + nnz(keep);
            end
        end

        keepSelf = ~isZeroShift;
        Rself = Rshifts(keepSelf, :);
        if ~isempty(Rself)
            x2 = sum(Rself.^2, 2);
            keep = (x2 > 0) & (x2 <= rcut2);
            nSelfPerSite = nnz(keep);
            nKeep = nKeep + nTarget * nSelfPerSite;
        end

    else
        nTarget = numel(targetSites);
        nSource = numel(sourceSites);

        for aa = 1:nTarget
            i = targetSites(aa);
            ri = pos(i, :);

            for bb = 1:nSource
                j = sourceSites(bb);
                rij0 = pos(j, :) - ri;

                xvec = rij0 + Rshifts;
                x2 = sum(xvec.^2, 2);

                if i == j
                    keep = (x2 > 0) & (x2 <= rcut2);
                else
                    keep = (x2 <= rcut2);
                end
                nKeep = nKeep + nnz(keep);
            end
        end
    end

    % Allocate exactly
    pair_i = zeros(nKeep, 1);
    pair_j = zeros(nKeep, 1);
    image_n = zeros(nKeep, 3);
    dr_all = zeros(nKeep, 3);
    r_bare_all = zeros(nKeep, 1);
    r2_bare_all = zeros(nKeep, 1);

    % PASS 2: fill arrays
    nFill = 0;

    if isSymmetric
        nTarget = numel(targetSites);

        for aa = 1:(nTarget - 1)
            i = targetSites(aa);
            ri = pos(i, :);

            for bb = (aa + 1):nTarget
                j = targetSites(bb);
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

                m = numel(xkeep);
                idx = (nFill + 1):(nFill + m);

                pair_i(idx) = i;
                pair_j(idx) = j;
                image_n(idx, :) = imageKeep;
                dr_all(idx, :) = xveckeep;
                r_bare_all(idx) = xkeep;
                r2_bare_all(idx) = x2keep;

                nFill = nFill + m;
            end
        end

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

                for aa = 1:nTarget
                    i = targetSites(aa);

                    m = numel(x);
                    idx = (nFill + 1):(nFill + m);

                    pair_i(idx) = i;
                    pair_j(idx) = i;
                    image_n(idx, :) = imageSelf;
                    dr_all(idx, :) = Rself;
                    r_bare_all(idx) = x;
                    r2_bare_all(idx) = x2;

                    nFill = nFill + m;
                end
            end
        end

    else
        nTarget = numel(targetSites);
        nSource = numel(sourceSites);

        for aa = 1:nTarget
            i = targetSites(aa);
            ri = pos(i, :);

            for bb = 1:nSource
                j = sourceSites(bb);
                rij0 = pos(j, :) - ri;

                xvec = rij0 + Rshifts;
                x2 = sum(xvec.^2, 2);

                if i == j
                    keep = (x2 > 0) & (x2 <= rcut2);
                else
                    keep = (x2 <= rcut2);
                end

                if ~any(keep)
                    continue;
                end

                x2keep = x2(keep);
                xkeep = sqrt(x2keep);
                xveckeep = xvec(keep, :);
                imageKeep = imageN(keep, :);

                m = numel(xkeep);
                idx = (nFill + 1):(nFill + m);

                pair_i(idx) = i;
                pair_j(idx) = j;
                image_n(idx, :) = imageKeep;
                dr_all(idx, :) = xveckeep;
                r_bare_all(idx) = xkeep;
                r2_bare_all(idx) = x2keep;

                nFill = nFill + m;
            end
        end
    end

    if nFill ~= nKeep
        error('geom:build_periodic_realspace_geom_cache:FillCountMismatch', ...
            'Two-pass fill count mismatch: counted %d, filled %d.', nKeep, nFill);
    end

    cache = struct();
    cache.mode = 'periodic_realspace_geom';
    cache.nSites = nSites;

    cache.target_mask = targetMask;
    cache.source_mask = sourceMask;
    cache.target_idx = targetSites;
    cache.source_idx = sourceSites;

    % Backward compatibility
    cache.site_mask = targetMask;
    cache.site_idx = targetSites;
    cache.activeSites = targetSites;

    cache.pair_i = pair_i;
    cache.pair_j = pair_j;
    cache.image_n = image_n;
    cache.dr = dr_all;

    cache.r_bare = r_bare_all;
    cache.r2_bare = r2_bare_all;

    cache.rcut = rcut;
    cache.real_image_bounds = [nxmax, nymax, nzmax];
    cache.nImageShifts = nR;
    cache.nInteractions = nKeep;
    cache.includes_self_image = true;
    cache.is_symmetric_pair_cache = isSymmetric;
end

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('geom:build_periodic_realspace_geom_cache:MissingLattice', ...
            'sys.super_lattice or sys.lattice is required.');
    end

    validateattributes(H, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
        mfilename, 'lattice');
end