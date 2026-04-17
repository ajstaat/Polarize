function rowCache = build_active_row_cache(sys, problem, pairCache)
%BUILD_ACTIVE_ROW_CACHE Build row-wise active-space neighbor cache for matrix-free GS/SOR.
%
% rowCache = geom.build_active_row_cache(sys, problem, pairCache)
%
% Inputs
%   sys       canonical polarization-system struct
%   problem   struct from thole.prepare_scf_problem(...)
%   pairCache unordered full-space pair cache from
%             geom.build_nonperiodic_pair_cache(...)
%
% Output
%   rowCache  struct with active-space row-wise neighbor data
%
% Fields
%   .nPolSites
%   .activeSites        full-site indices for active polarizable sites
%   .full_to_active     map full site index -> active row index (0 if inactive)
%   .row_ptr            nPolSites+1 x 1 CSR row pointer
%   .col_idx            directed active-space neighbor columns
%   .dr                 directed displacement vectors r_j - r_i
%   .r_bare
%   .r2_bare
%   .inv_r3_bare
%   .inv_r5_bare
%   .thole_f3
%   .thole_f5
%
% Notes
%   - Each unordered pair in pairCache generates two directed row entries:
%       i <- j   and   j <- i
%   - dr is always stored as r_j - r_i for the directed interaction i <- j.
%   - This layout is intended for row-wise GS/SOR sweeps.

    narginchk(3, 3);

    if ~isfield(problem, 'activeSites') || isempty(problem.activeSites)
        error('geom:build_active_row_cache:MissingActiveSites', ...
            'problem.activeSites is required.');
    end

    activeSites = problem.activeSites(:);
    nPolSites = numel(activeSites);
    nSites = problem.nSites;

    if pairCache.nSites ~= nSites
        error('geom:build_active_row_cache:SiteCountMismatch', ...
            'pairCache.nSites must match problem.nSites.');
    end

    fullToActive = zeros(nSites, 1);
    fullToActive(activeSites) = 1:nPolSites;

    pair_i_full = pairCache.pair_i(:);
    pair_j_full = pairCache.pair_j(:);

    ai = fullToActive(pair_i_full);
    aj = fullToActive(pair_j_full);

    keep = (ai > 0) & (aj > 0);
    ai = ai(keep);
    aj = aj(keep);

    dr_pair      = pairCache.dr(keep, :);
    r_bare_pair  = pairCache.r_bare(keep);
    r2_bare_pair = pairCache.r2_bare(keep);
    invr3_pair   = pairCache.inv_r3_bare(keep);
    invr5_pair   = pairCache.inv_r5_bare(keep);

    haveThole = isfield(pairCache, 'thole_f3') && isfield(pairCache, 'thole_f5');
    if haveThole
        f3_pair = pairCache.thole_f3(keep);
        f5_pair = pairCache.thole_f5(keep);
    else
        f3_pair = [];
        f5_pair = [];
    end

    nPairs = numel(ai);
    nDirected = 2 * nPairs;

    row_idx  = zeros(nDirected, 1);
    col_idx  = zeros(nDirected, 1);
    dr       = zeros(nDirected, 3);
    r_bare   = zeros(nDirected, 1);
    r2_bare  = zeros(nDirected, 1);
    inv_r3   = zeros(nDirected, 1);
    inv_r5   = zeros(nDirected, 1);

    if haveThole
        thole_f3 = zeros(nDirected, 1);
        thole_f5 = zeros(nDirected, 1);
    else
        thole_f3 = [];
        thole_f5 = [];
    end

    idx = 0;
    for p = 1:nPairs
        % Directed entry: ai <- aj
        idx = idx + 1;
        row_idx(idx) = ai(p);
        col_idx(idx) = aj(p);
        dr(idx, :) = dr_pair(p, :);        % r_j - r_i
        r_bare(idx) = r_bare_pair(p);
        r2_bare(idx) = r2_bare_pair(p);
        inv_r3(idx) = invr3_pair(p);
        inv_r5(idx) = invr5_pair(p);
        if haveThole
            thole_f3(idx) = f3_pair(p);
            thole_f5(idx) = f5_pair(p);
        end

        % Directed entry: aj <- ai
        idx = idx + 1;
        row_idx(idx) = aj(p);
        col_idx(idx) = ai(p);
        dr(idx, :) = -dr_pair(p, :);       % now r_i - r_j for target j <- source i
        r_bare(idx) = r_bare_pair(p);
        r2_bare(idx) = r2_bare_pair(p);
        inv_r3(idx) = invr3_pair(p);
        inv_r5(idx) = invr5_pair(p);
        if haveThole
            thole_f3(idx) = f3_pair(p);
            thole_f5(idx) = f5_pair(p);
        end
    end

    % Sort by row for CSR layout
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