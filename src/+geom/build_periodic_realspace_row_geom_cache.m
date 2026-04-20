function rowCache = build_periodic_realspace_row_geom_cache(sys, problem, geomCache)
%BUILD_PERIODIC_REALSPACE_ROW_GEOM_CACHE Build directed row cache
% from periodic real-space geometry cache.
%
% rowCache = geom.build_periodic_realspace_row_geom_cache(sys, problem, geomCache)
%
% Supports both:
%   - symmetric unordered dipole-style geometry caches
%   - directed rectangular target/source charge-style geometry caches
%
% Output fields
%   .mode
%   .nSites
%   .nTargetSites
%   .nSourceSites
%   .targetSites
%   .sourceSites
%   .activeSites        alias of targetSites for backward compatibility
%   .full_to_target
%   .full_to_active     alias of full_to_target
%   .full_to_source
%   .row_ptr
%   .col_idx            source-local index
%   .source_full_idx
%   .dr
%   .r_bare
%   .r2_bare
%   .rcut
%   .nInteractions

    narginchk(3, 3);

    validateattributes(problem, {'struct'}, {'scalar'}, mfilename, 'problem', 2);
    validateattributes(geomCache, {'struct'}, {'scalar'}, mfilename, 'geomCache', 3);

    nSites = problem.nSites;
    if ~isfield(geomCache, 'nSites') || geomCache.nSites ~= nSites
        error('geom:build_periodic_realspace_row_geom_cache:SiteCountMismatch', ...
            'geomCache.nSites must match problem.nSites.');
    end

    if ~isfield(geomCache, 'target_idx') || isempty(geomCache.target_idx)
        error('geom:build_periodic_realspace_row_geom_cache:MissingTargetIdx', ...
            'geomCache.target_idx is required.');
    end
    if ~isfield(geomCache, 'source_idx') || isempty(geomCache.source_idx)
        error('geom:build_periodic_realspace_row_geom_cache:MissingSourceIdx', ...
            'geomCache.source_idx is required.');
    end

    targetSites = geomCache.target_idx(:);
    sourceSites = geomCache.source_idx(:);

    nTargetSites = numel(targetSites);
    nSourceSites = numel(sourceSites);

    fullToTarget = zeros(nSites, 1);
    fullToTarget(targetSites) = 1:nTargetSites;

    fullToSource = zeros(nSites, 1);
    fullToSource(sourceSites) = 1:nSourceSites;

    pair_i_full = geomCache.pair_i(:);
    pair_j_full = geomCache.pair_j(:);
    dr_pair = geomCache.dr;
    r_bare_pair = geomCache.r_bare(:);
    r2_bare_pair = geomCache.r2_bare(:);

    ti = fullToTarget(pair_i_full);
    sj = fullToSource(pair_j_full);

    keep = (ti > 0) & (sj > 0);
    if ~any(keep)
        rowCache = struct();
        rowCache.mode = 'periodic_realspace_row_geom';
        rowCache.nSites = nSites;
        rowCache.nTargetSites = nTargetSites;
        rowCache.nSourceSites = nSourceSites;
        rowCache.nPolSites = nTargetSites; % backward compatibility
        rowCache.targetSites = targetSites;
        rowCache.sourceSites = sourceSites;
        rowCache.activeSites = targetSites;
        rowCache.full_to_target = fullToTarget;
        rowCache.full_to_active = fullToTarget;
        rowCache.full_to_source = fullToSource;
        rowCache.row_ptr = ones(nTargetSites + 1, 1);
        rowCache.col_idx = zeros(0, 1);
        rowCache.source_full_idx = zeros(0, 1);
        rowCache.dr = zeros(0, 3);
        rowCache.r_bare = zeros(0, 1);
        rowCache.r2_bare = zeros(0, 1);
        rowCache.rcut = geomCache.rcut;
        rowCache.nInteractions = 0;
        return;
    end

    pair_i_full = pair_i_full(keep);
    pair_j_full = pair_j_full(keep);
    ti = ti(keep);
    sj = sj(keep);

    dr_pair = dr_pair(keep, :);
    r_bare_pair = r_bare_pair(keep);
    r2_bare_pair = r2_bare_pair(keep);

    isSymmetric = isfield(geomCache, 'is_symmetric_pair_cache') && geomCache.is_symmetric_pair_cache;

    if isSymmetric
        % Expand unordered cache into directed row entries.
        nEntries = numel(ti);
        isSelf = (pair_i_full == pair_j_full);
        nDirected = 2 * nnz(~isSelf) + nnz(isSelf);

        row_idx = zeros(nDirected, 1);
        col_idx = zeros(nDirected, 1);
        source_full_idx = zeros(nDirected, 1);
        dr = zeros(nDirected, 3);
        r_bare = zeros(nDirected, 1);
        r2_bare = zeros(nDirected, 1);

        idx = 0;
        for p = 1:nEntries
            if isSelf(p)
                idx = idx + 1;
                row_idx(idx) = ti(p);
                col_idx(idx) = sj(p);
                source_full_idx(idx) = pair_i_full(p);

                dr(idx, :) = dr_pair(p, :);
                r_bare(idx) = r_bare_pair(p);
                r2_bare(idx) = r2_bare_pair(p);
            else
                % i <- j
                idx = idx + 1;
                row_idx(idx) = ti(p);
                col_idx(idx) = sj(p);
                source_full_idx(idx) = pair_j_full(p);

                dr(idx, :) = dr_pair(p, :);
                r_bare(idx) = r_bare_pair(p);
                r2_bare(idx) = r2_bare_pair(p);

                % j <- i
                idx = idx + 1;
                row_idx(idx) = fullToTarget(pair_j_full(p));
                col_idx(idx) = fullToSource(pair_i_full(p));
                source_full_idx(idx) = pair_i_full(p);

                dr(idx, :) = -dr_pair(p, :);
                r_bare(idx) = r_bare_pair(p);
                r2_bare(idx) = r2_bare_pair(p);
            end
        end

    else
        % Already directed target <- source geometry.
        nDirected = numel(ti);

        row_idx = ti;
        col_idx = sj;
        source_full_idx = pair_j_full;
        dr = dr_pair;
        r_bare = r_bare_pair;
        r2_bare = r2_bare_pair;
    end

    [row_idx, perm] = sort(row_idx, 'ascend');
    col_idx = col_idx(perm);
    source_full_idx = source_full_idx(perm);
    dr = dr(perm, :);
    r_bare = r_bare(perm);
    r2_bare = r2_bare(perm);

    row_ptr = zeros(nTargetSites + 1, 1);
    cursor = 1;

    for i = 1:nTargetSites
        row_ptr(i) = cursor;
        while cursor <= nDirected && row_idx(cursor) == i
            cursor = cursor + 1;
        end
        row_ptr(i + 1) = cursor;
    end

    rowCache = struct();
    rowCache.mode = 'periodic_realspace_row_geom';
    rowCache.nSites = nSites;
    rowCache.nTargetSites = nTargetSites;
    rowCache.nSourceSites = nSourceSites;
    rowCache.nPolSites = nTargetSites; % backward compatibility

    rowCache.targetSites = targetSites;
    rowCache.sourceSites = sourceSites;
    rowCache.activeSites = targetSites; % backward compatibility

    rowCache.full_to_target = fullToTarget;
    rowCache.full_to_active = fullToTarget; % backward compatibility
    rowCache.full_to_source = fullToSource;

    rowCache.row_ptr = row_ptr;
    rowCache.col_idx = col_idx;
    rowCache.source_full_idx = source_full_idx;

    rowCache.dr = dr;
    rowCache.r_bare = r_bare;
    rowCache.r2_bare = r2_bare;

    rowCache.rcut = geomCache.rcut;
    rowCache.nInteractions = nDirected;
end