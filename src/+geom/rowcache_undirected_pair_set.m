function [P, counts, rows, cols] = rowcache_undirected_pair_set(rowCache)
%ROWCACHE_UNDIRECTED_PAIR_SET Unique undirected pair set from a row cache.
%
%   P = geom.rowcache_undirected_pair_set(rowCache)
%   [P, counts] = geom.rowcache_undirected_pair_set(rowCache)
%   [P, counts, rows, cols] = geom.rowcache_undirected_pair_set(rowCache)
%
% Inputs
%   rowCache : struct with fields
%              .row_ptr
%              .col_idx
%              .nActive
%
% Outputs
%   P      : nPairsUnique x 2 array of sorted undirected pairs [i j], i < j
%   counts : nPairsUnique x 1 multiplicity of each undirected pair in the
%            directed row cache, after collapsing (i,j) and (j,i)
%   rows   : nDir x 1 source-row index reconstructed from CSR row_ptr
%   cols   : nDir x 1 copied from rowCache.col_idx
%
% Notes
%   - Reconstructs the directed row->col edge list from CSR storage.
%   - Converts each directed edge to undirected form [min(i,j) max(i,j)].
%   - Removes self edges.
%   - Returns unique undirected pairs in deterministic sorted order.
%   - Optional counts are useful for diagnosing duplicate insertion.

if ~isfield(rowCache, 'row_ptr') || isempty(rowCache.row_ptr)
    error('geom:rowcache_undirected_pair_set:MissingRowPtr', ...
        'rowCache.row_ptr is required.');
end
if ~isfield(rowCache, 'col_idx') || isempty(rowCache.col_idx)
    error('geom:rowcache_undirected_pair_set:MissingColIdx', ...
        'rowCache.col_idx is required.');
end
if ~isfield(rowCache, 'nActive') || isempty(rowCache.nActive)
    error('geom:rowcache_undirected_pair_set:MissingNActive', ...
        'rowCache.nActive is required.');
end

nActive = rowCache.nActive;
row_ptr = rowCache.row_ptr(:);
cols = rowCache.col_idx(:);

if numel(row_ptr) ~= nActive + 1
    error('geom:rowcache_undirected_pair_set:BadRowPtrSize', ...
        'row_ptr must have length nActive + 1.');
end

nDir = numel(cols);
rows = zeros(nDir, 1);

for i = 1:nActive
    k0 = row_ptr(i);
    k1 = row_ptr(i+1) - 1;
    if k1 >= k0
        rows(k0:k1) = i;
    end
end

% Build undirected edge list
E = [min(rows, cols), max(rows, cols)];

% Drop diagonal/self edges if any slipped in
mask = E(:,1) < E(:,2);
E = E(mask, :);

if isempty(E)
    P = zeros(0, 2);
    counts = zeros(0, 1);
    return;
end

% Deterministic sorted unique pairs with multiplicities
[P, ~, ic] = unique(E, 'rows', 'sorted');
counts = accumarray(ic, 1, [size(P,1), 1], @sum, 0);
end