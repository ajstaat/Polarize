function pairs = query_pairs_within_cutoff(spatial, cutoff, opts)
%QUERY_PAIRS_WITHIN_CUTOFF Find unordered pairs within a distance cutoff.
%
% pairs = geom.query_pairs_within_cutoff(spatial, cutoff)
% pairs = geom.query_pairs_within_cutoff(spatial, cutoff, opts)
%
% Inputs
%   spatial : struct from geom.build_spatial_index(...)
%   cutoff  : positive scalar cutoff distance
%   opts    : optional struct with fields
%             .subset_idx vector of point indices into spatial.pos
%                        default = all points
%             .return_r logical, default true
%             .return_dr logical, default true
%             .return_full_idx logical, default false
%
% Output
%   pairs : struct with fields
%           .i, .j unordered pair indices
%           .r2 squared distances
%           .dr displacement vectors, if requested
%           .r distances, if requested
%           .index_space 'subset' or 'full'
%           .subset_idx subset mapping used for the query
%
% Notes
% - If subset_idx is supplied and return_full_idx = false, i and j are in
%   subset-local indexing (1:numel(subset_idx)).
% - If subset_idx is supplied and return_full_idx = true, i and j are in
%   full indexing into spatial.pos.
% - For periodic geometry, minimum-image displacement is used.

narginchk(2, 3);
if nargin < 3 || isempty(opts)
    opts = struct();
end

validateattributes(cutoff, {'double'}, {'scalar', 'positive', 'finite', 'real'}, ...
    mfilename, 'cutoff', 2);

if ~isfield(spatial, 'method') || isempty(spatial.method)
    error('geom:query_pairs_within_cutoff:BadSpatial', ...
        'spatial.method is missing.');
end

switch char(string(spatial.method))
    case 'bruteforce'
        pairs = local_query_pairs_bruteforce(spatial, cutoff, opts);
    case 'cell_list'
        pairs = geom.query_pairs_cell_list(spatial, cutoff, opts);
    otherwise
        error('geom:query_pairs_within_cutoff:UnsupportedMethod', ...
            'Unsupported spatial.method "%s".', spatial.method);
end

end

function pairs = local_query_pairs_bruteforce(spatial, cutoff, opts)

subsetIdx     = local_get_opt(opts, 'subset_idx', []);
returnR       = local_get_opt(opts, 'return_r', true);
returnDr      = local_get_opt(opts, 'return_dr', true);
returnFullIdx = local_get_opt(opts, 'return_full_idx', false);

if isempty(subsetIdx)
    idxMap = (1:spatial.n).';
else
    validateattributes(subsetIdx, {'numeric'}, {'vector', 'integer', 'positive'}, ...
        mfilename, 'opts.subset_idx');
    idxMap = subsetIdx(:);
    if any(idxMap > spatial.n)
        error('geom:query_pairs_within_cutoff:SubsetOutOfRange', ...
            'opts.subset_idx contains indices outside spatial.pos.');
    end
end

pos = spatial.pos(idxMap, :);
n = size(pos, 1);
cutoff2 = cutoff^2;

maxPairs = n * (n - 1) / 2;
pairI = zeros(maxPairs, 1);
pairJ = zeros(maxPairs, 1);
r2All = zeros(maxPairs, 1);
if returnDr
    drAll = zeros(maxPairs, 3);
else
    drAll = zeros(0, 3);
end

nPairs = 0;
for i = 1:(n-1)
    ri = pos(i, :);
    for j = (i+1):n
        dr = pos(j, :) - ri;
        if spatial.isPeriodic
            dr = local_minimum_image(dr, spatial.cell, spatial.invCell);
        end
        r2 = dot(dr, dr);
        if r2 <= cutoff2
            nPairs = nPairs + 1;
            pairI(nPairs) = i;
            pairJ(nPairs) = j;
            r2All(nPairs) = r2;
            if returnDr
                drAll(nPairs, :) = dr;
            end
        end
    end
end

pairI = pairI(1:nPairs);
pairJ = pairJ(1:nPairs);
r2All = r2All(1:nPairs);
if returnDr
    drAll = drAll(1:nPairs, :);
end

if returnFullIdx
    pairI = idxMap(pairI);
    pairJ = idxMap(pairJ);
    indexSpace = 'full';
else
    indexSpace = 'subset';
end

pairs = struct();
pairs.i = pairI;
pairs.j = pairJ;
pairs.r2 = r2All;
if returnDr
    pairs.dr = drAll;
end
if returnR
    pairs.r = sqrt(r2All);
end
pairs.index_space = indexSpace;
pairs.subset_idx = idxMap;

end

function dr = local_minimum_image(dr, cellMat, invCell)
frac = (invCell * dr(:)).';
frac = frac - round(frac);
dr = (cellMat * frac(:)).';
end

function value = local_get_opt(s, name, defaultValue)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end