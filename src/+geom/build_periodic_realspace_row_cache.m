function rowCache = build_periodic_realspace_row_cache(sys, problem, realCache)
%BUILD_PERIODIC_REALSPACE_ROW_CACHE Build directed active-space row cache
% for periodic real-space Ewald interactions.
%
% rowCache = ewald.build_periodic_realspace_row_cache(sys, problem, realCache)
%
% Inputs
%   sys        canonical polarization-system struct
%   problem    struct from thole.prepare_scf_problem(...)
%   realCache  unordered periodic real-space cache from
%              ewald.build_periodic_realspace_cache(...)
%
% Output
%   rowCache   struct with active-space CSR-style directed row data
%
% Fields
%   .mode               'periodic_realspace_row'
%   .nSites
%   .nPolSites
%   .activeSites
%   .full_to_active
%   .row_ptr            (nPolSites+1) x 1 CSR row pointer
%   .col_idx            directed active-space source columns
%   .source_full_idx    directed full-site source indices
%   .dr                 directed displacement vectors for target i <- source j
%   .r_bare
%   .r2_bare
%   .coeff_iso
%   .coeff_dyad
%   .alpha
%   .rcut
%   .use_thole
%   .nInteractions
%
% Notes
%   - Distinct-site entries in realCache are stored once as unordered
%     interactions. Each such entry is expanded here into two directed row
%     entries:
%         i <- j   with dr = r_j + R - r_i
%         j <- i   with dr = -(r_j + R - r_i)
%   - Self-image entries (i == j, R ~= 0) remain single directed row
%     entries on their own row.
%   - The tensor action used later is:
%         E_{i<-j} = coeff_iso * mu_j + coeff_dyad * dr * (dr · mu_j)

    narginchk(3, 3);

    validateattributes(problem, {'struct'}, {'scalar'}, mfilename, 'problem', 2);
    validateattributes(realCache, {'struct'}, {'scalar'}, mfilename, 'realCache', 3);

    if ~isfield(problem, 'activeSites') || isempty(problem.activeSites)
        error('ewald:build_periodic_realspace_row_cache:MissingActiveSites', ...
            'problem.activeSites is required.');
    end

    activeSites = problem.activeSites(:);
    nPolSites = numel(activeSites);
    nSites = problem.nSites;

    if ~isfield(realCache, 'nSites') || realCache.nSites ~= nSites
        error('ewald:build_periodic_realspace_row_cache:SiteCountMismatch', ...
            'realCache.nSites must match problem.nSites.');
    end

    if ~isfield(realCache, 'pair_i') || ~isfield(realCache, 'pair_j') || ...
       ~isfield(realCache, 'dr') || ~isfield(realCache, 'coeff_iso') || ...
       ~isfield(realCache, 'coeff_dyad')
        error('ewald:build_periodic_realspace_row_cache:MissingCacheFields', ...
            'realCache is missing required interaction fields.');
    end

    fullToActive = zeros(nSites, 1);
    fullToActive(activeSites) = 1:nPolSites;

    pair_i_full = realCache.pair_i(:);
    pair_j_full = realCache.pair_j(:);

    ai = fullToActive(pair_i_full);
    aj = fullToActive(pair_j_full);

    keep = (ai > 0) & (aj > 0);
    if ~any(keep)
        % Return an empty-but-valid row cache.
        rowCache = struct();
        rowCache.mode = 'periodic_realspace_row';
        rowCache.nSites = nSites;
        rowCache.nPolSites = nPolSites;
        rowCache.activeSites = activeSites;
        rowCache.full_to_active = fullToActive;
        rowCache.row_ptr = ones(nPolSites + 1, 1);
        rowCache.col_idx = zeros(0, 1);
        rowCache.source_full_idx = zeros(0, 1);
        rowCache.dr = zeros(0, 3);
        rowCache.r_bare = zeros(0, 1);
        rowCache.r2_bare = zeros(0, 1);
        rowCache.coeff_iso = zeros(0, 1);
        rowCache.coeff_dyad = zeros(0, 1);
        rowCache.alpha = realCache.alpha;
        rowCache.rcut = realCache.rcut;
        rowCache.use_thole = realCache.use_thole;
        rowCache.nInteractions = 0;
        return;
    end

    pair_i_full = pair_i_full(keep);
    pair_j_full = pair_j_full(keep);
    ai = ai(keep);
    aj = aj(keep);

    dr_pair      = realCache.dr(keep, :);
    r_bare_pair  = realCache.r_bare(keep);
    r2_bare_pair = realCache.r2_bare(keep);
    coeff_iso_pair  = realCache.coeff_iso(keep);
    coeff_dyad_pair = realCache.coeff_dyad(keep);

    nEntries = numel(ai);

    % Count directed entries:
    %   distinct-site -> two directed entries
    %   self-image    -> one directed entry
    isSelf = (pair_i_full == pair_j_full);

    nDirected = 2 * nnz(~isSelf) + nnz(isSelf);

    row_idx = zeros(nDirected, 1);
    col_idx = zeros(nDirected, 1);
    source_full_idx = zeros(nDirected, 1);
    dr = zeros(nDirected, 3);
    r_bare = zeros(nDirected, 1);
    r2_bare = zeros(nDirected, 1);
    coeff_iso = zeros(nDirected, 1);
    coeff_dyad = zeros(nDirected, 1);

    idx = 0;
    for p = 1:nEntries
        if isSelf(p)
            % Self-image entry stays on one directed row.
            idx = idx + 1;
            row_idx(idx) = ai(p);
            col_idx(idx) = ai(p);
            source_full_idx(idx) = pair_i_full(p);

            dr(idx, :) = dr_pair(p, :);
            r_bare(idx) = r_bare_pair(p);
            r2_bare(idx) = r2_bare_pair(p);
            coeff_iso(idx) = coeff_iso_pair(p);
            coeff_dyad(idx) = coeff_dyad_pair(p);
        else
            % Directed entry: ai <- aj
            idx = idx + 1;
            row_idx(idx) = ai(p);
            col_idx(idx) = aj(p);
            source_full_idx(idx) = pair_j_full(p);

            dr(idx, :) = dr_pair(p, :);      % r_j + R - r_i
            r_bare(idx) = r_bare_pair(p);
            r2_bare(idx) = r2_bare_pair(p);
            coeff_iso(idx) = coeff_iso_pair(p);
            coeff_dyad(idx) = coeff_dyad_pair(p);

            % Directed entry: aj <- ai
            idx = idx + 1;
            row_idx(idx) = aj(p);
            col_idx(idx) = ai(p);
            source_full_idx(idx) = pair_i_full(p);

            dr(idx, :) = -dr_pair(p, :);     % r_i - R - r_j = -(r_j + R - r_i)
            r_bare(idx) = r_bare_pair(p);
            r2_bare(idx) = r2_bare_pair(p);
            coeff_iso(idx) = coeff_iso_pair(p);
            coeff_dyad(idx) = coeff_dyad_pair(p);
        end
    end

    % Sort by row for CSR layout.
    [row_idx, perm] = sort(row_idx, 'ascend');
    col_idx = col_idx(perm);
    source_full_idx = source_full_idx(perm);
    dr = dr(perm, :);
    r_bare = r_bare(perm);
    r2_bare = r2_bare(perm);
    coeff_iso = coeff_iso(perm);
    coeff_dyad = coeff_dyad(perm);

    row_ptr = zeros(nPolSites + 1, 1);
    cursor = 1;

    for i = 1:nPolSites
        row_ptr(i) = cursor;
        while cursor <= nDirected && row_idx(cursor) == i
            cursor = cursor + 1;
        end
        row_ptr(i + 1) = cursor;
    end

    rowCache = struct();
    rowCache.mode = 'periodic_realspace_row';
    rowCache.nSites = nSites;
    rowCache.nPolSites = nPolSites;
    rowCache.activeSites = activeSites;
    rowCache.full_to_active = fullToActive;

    rowCache.row_ptr = row_ptr;
    rowCache.col_idx = col_idx;
    rowCache.source_full_idx = source_full_idx;

    rowCache.dr = dr;
    rowCache.r_bare = r_bare;
    rowCache.r2_bare = r2_bare;

    rowCache.coeff_iso = coeff_iso;
    rowCache.coeff_dyad = coeff_dyad;

    rowCache.alpha = realCache.alpha;
    rowCache.rcut = realCache.rcut;
    rowCache.use_thole = realCache.use_thole;
    rowCache.nInteractions = nDirected;
end