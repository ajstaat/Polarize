function rowCache = build_active_row_cache_direct(sys, problem, opts, spatial)
%BUILD_ACTIVE_ROW_CACHE_DIRECT Build row-wise active-space neighbor cache
%directly, without first materializing an unordered full-space pair cache.
%
% rowCache = geom.build_active_row_cache_direct(sys, problem)
% rowCache = geom.build_active_row_cache_direct(sys, problem, opts)
% rowCache = geom.build_active_row_cache_direct(sys, problem, opts, spatial)
%
% Inputs
%   sys       canonical polarization-system struct
%   problem   struct from thole.prepare_scf_problem(...)
%   opts      optional struct with fields:
%       .rcut    positive scalar cutoff in bohr, default = Inf
%
%   spatial   optional spatial index built on ACTIVE-SITE positions
%
% Output
%   rowCache  struct with active-space row-wise neighbor data
%
% Fields
%   .nPolSites
%   .activeSites
%   .full_to_active
%   .row_ptr
%   .col_idx
%   .dr
%   .r_bare
%   .r2_bare
%   .inv_r3_bare
%   .inv_r5_bare
%   .thole_f3
%   .thole_f5
%
% Notes
%   - Builds the final directed active-space row structure used by
%     matrix-free GS/SOR directly.
%   - dr is stored as r_j - r_i for directed interaction i <- j.

    narginchk(2, 4);

    if nargin < 3 || isempty(opts)
        opts = struct();
    end

    if ~isfield(problem, 'activeSites') || isempty(problem.activeSites)
        error('geom:build_active_row_cache_direct:MissingActiveSites', ...
            'problem.activeSites is required.');
    end

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('geom:build_active_row_cache_direct:MissingSitePos', ...
            'sys.site_pos is required.');
    end

    activeSites = problem.activeSites(:);
    nPolSites = numel(activeSites);
    nSites = problem.nSites;

    fullToActive = zeros(nSites, 1);
    fullToActive(activeSites) = 1:nPolSites;

    rcut = inf;
    if isfield(opts, 'rcut') && ~isempty(opts.rcut)
        rcut = opts.rcut;
    end
    if ~(isscalar(rcut) && isreal(rcut) && isfinite(rcut) && rcut > 0) && ~isinf(rcut)
        error('geom:build_active_row_cache_direct:BadRcut', ...
            'opts.rcut must be a positive scalar or Inf.');
    end

    posAct = sys.site_pos(activeSites, :);

    if nargin < 4 || isempty(spatial)
        spatial = geom.build_spatial_index(posAct, struct('isPeriodic', false));
    end

    queryOpts = struct();
    queryOpts.subset_idx = (1:nPolSites).';
    queryOpts.return_r = true;
    queryOpts.return_dr = true;
    queryOpts.return_full_idx = true;

    pairs = geom.query_pairs_within_cutoff(spatial, rcut, queryOpts);

    pair_i_local = pairs.i(:);
    pair_j_local = pairs.j(:);
    dr_pair      = pairs.dr;
    r_bare_pair  = pairs.r(:);
    r2_bare_pair = pairs.r2(:);

    if any(r2_bare_pair <= 0)
        error('geom:build_active_row_cache_direct:NonPositiveDistance', ...
            'Encountered non-positive pair distance while building row cache.');
    end

    invr3_pair = 1 ./ (r_bare_pair.^3);
    invr5_pair = invr3_pair ./ r2_bare_pair;

    haveThole = isfield(sys, 'site_alpha') && ~isempty(sys.site_alpha) && ...
                isfield(sys, 'thole_a')   && ~isempty(sys.thole_a);

    if haveThole
        alpha_full = sys.site_alpha(:);
        if numel(alpha_full) ~= nSites
            error('geom:build_active_row_cache_direct:BadSiteAlpha', ...
                'sys.site_alpha must have length problem.nSites.');
        end
        alpha_act = alpha_full(activeSites);
        a = sys.thole_a;

        nPairs = numel(pair_i_local);
        f3_pair = zeros(nPairs, 1);
        f5_pair = zeros(nPairs, 1);
        for p = 1:nPairs
            tf = thole.thole_f3f5_factors(r_bare_pair(p), ...
                alpha_act(pair_i_local(p)), alpha_act(pair_j_local(p)), a);
            f3_pair(p) = tf.f3;
            f5_pair(p) = tf.f5;
        end
    else
        f3_pair = [];
        f5_pair = [];
    end

    nPairs = numel(pair_i_local);
    nDirected = 2 * nPairs;

    row_idx = zeros(nDirected, 1);
    col_idx = zeros(nDirected, 1);
    dr = zeros(nDirected, 3);
    r_bare = zeros(nDirected, 1);
    r2_bare = zeros(nDirected, 1);
    inv_r3 = zeros(nDirected, 1);
    inv_r5 = zeros(nDirected, 1);

    if haveThole
        thole_f3 = zeros(nDirected, 1);
        thole_f5 = zeros(nDirected, 1);
    else
        thole_f3 = [];
        thole_f5 = [];
    end

    idx = 0;
    for p = 1:nPairs
        % Directed entry: i <- j
        idx = idx + 1;
        row_idx(idx) = pair_i_local(p);
        col_idx(idx) = pair_j_local(p);
        dr(idx, :) = dr_pair(p, :);
        r_bare(idx) = r_bare_pair(p);
        r2_bare(idx) = r2_bare_pair(p);
        inv_r3(idx) = invr3_pair(p);
        inv_r5(idx) = invr5_pair(p);
        if haveThole
            thole_f3(idx) = f3_pair(p);
            thole_f5(idx) = f5_pair(p);
        end

        % Directed entry: j <- i
        idx = idx + 1;
        row_idx(idx) = pair_j_local(p);
        col_idx(idx) = pair_i_local(p);
        dr(idx, :) = -dr_pair(p, :);
        r_bare(idx) = r_bare_pair(p);
        r2_bare(idx) = r2_bare_pair(p);
        inv_r3(idx) = invr3_pair(p);
        inv_r5(idx) = invr5_pair(p);
        if haveThole
            thole_f3(idx) = f3_pair(p);
            thole_f5(idx) = f5_pair(p);
        end
    end

    [row_idx, perm] = sort(row_idx);
    col_idx = col_idx(perm);
    dr = dr(perm, :);
    r_bare = r_bare(perm);
    r2_bare = r2_bare(perm);
    inv_r3 = inv_r3(perm);
    inv_r5 = inv_r5(perm);
    if haveThole
        thole_f3 = thole_f3(perm);
        thole_f5 = thole_f5(perm);
    end

    row_ptr = zeros(nPolSites + 1, 1);
    row_ptr(1) = 1;

    cursor = 1;
    for i = 1:nPolSites
        while cursor <= nDirected && row_idx(cursor) < i
            cursor = cursor + 1;
        end
        row_ptr(i) = cursor;

        while cursor <= nDirected && row_idx(cursor) == i
            cursor = cursor + 1;
        end
        row_ptr(i + 1) = cursor;
    end

    rowCache = struct();
    rowCache.nPolSites = nPolSites;
    rowCache.activeSites = activeSites;
    rowCache.full_to_active = fullToActive;

    rowCache.row_ptr = row_ptr;
    rowCache.col_idx = col_idx;
    rowCache.dr = dr;

    rowCache.r_bare = r_bare;
    rowCache.r2_bare = r2_bare;
    rowCache.inv_r3_bare = inv_r3;
    rowCache.inv_r5_bare = inv_r5;

    if haveThole
        rowCache.thole_f3 = thole_f3;
        rowCache.thole_f5 = thole_f5;
    end
end