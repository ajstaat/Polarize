function pairs = query_pairs_cell_list(spatial, cutoff, opts)
%QUERY_PAIRS_CELL_LIST Query unordered pairs within cutoff from a cell list.
%
% pairs = geom.query_pairs_cell_list(spatial, cutoff)
% pairs = geom.query_pairs_cell_list(spatial, cutoff, opts)
%
% Inputs
%   spatial : struct from geom.build_cell_list(...) or
%             geom.build_spatial_index(..., method='cell_list')
%   cutoff  : positive scalar query cutoff, must satisfy cutoff <= spatial.cutoff
%   opts    : optional struct with fields
%             .subset_idx       indices into spatial.pos, default all
%             .return_r         logical, default true
%             .return_dr        logical, default true
%             .return_full_idx  logical, default false
%
% Output
%   pairs : struct with fields
%           .i, .j
%           .r2
%           .dr (if requested)
%           .r  (if requested)
%           .index_space
%           .subset_idx

narginchk(2, 3);
if nargin < 3 || isempty(opts)
    opts = struct();
end

validateattributes(cutoff, {'double'}, {'scalar', 'positive', 'finite', 'real'}, ...
    mfilename, 'cutoff', 2);

if cutoff > spatial.cutoff
    error('geom:query_pairs_cell_list:CutoffExceedsIndex', ...
        'Query cutoff exceeds spatial.cutoff used to build the cell list.');
end

subsetIdx     = local_get_opt(opts, 'subset_idx', []);
returnR       = local_get_opt(opts, 'return_r', true);
returnDr      = local_get_opt(opts, 'return_dr', true);
returnFullIdx = local_get_opt(opts, 'return_full_idx', false);

if isempty(subsetIdx)
    activeMask = true(spatial.n, 1);
    idxMap = (1:spatial.n).';
else
    validateattributes(subsetIdx, {'numeric'}, {'vector', 'integer', 'positive'}, ...
        mfilename, 'opts.subset_idx');
    idxMap = subsetIdx(:);
    if any(idxMap > spatial.n)
        error('geom:query_pairs_cell_list:SubsetOutOfRange', ...
            'opts.subset_idx contains indices outside spatial.pos.');
    end
    activeMask = false(spatial.n, 1);
    activeMask(idxMap) = true;
end

fullToLocal = zeros(spatial.n, 1);
fullToLocal(idxMap) = 1:numel(idxMap);

cutoff2 = cutoff^2;
gridShape = spatial.grid_shape;
offsets = spatial.neighbor_offsets;
nOffsets = size(offsets, 1);

alloc = max(1024, min(spatial.n * 32, max(1024, floor(spatial.n*(spatial.n-1)/2))));
pairI = zeros(alloc, 1);
pairJ = zeros(alloc, 1);
r2All = zeros(alloc, 1);
if returnDr
    drAll = zeros(alloc, 3);
else
    drAll = zeros(0, 3);
end
nPairs = 0;

for binLin = 1:spatial.n_bins
    iHead = spatial.bin_head(binLin);
    if iHead == 0
        continue;
    end

    [ix, iy, iz] = local_unlinearize(binLin, gridShape);

    for oo = 1:nOffsets
        off = offsets(oo, :);
        jx = ix + off(1);
        jy = iy + off(2);
        jz = iz + off(3);

        if spatial.isPeriodic
            jx = 1 + mod(jx - 1, gridShape(1));
            jy = 1 + mod(jy - 1, gridShape(2));
            jz = 1 + mod(jz - 1, gridShape(3));
        else
            if jx < 1 || jx > gridShape(1) || ...
               jy < 1 || jy > gridShape(2) || ...
               jz < 1 || jz > gridShape(3)
                continue;
            end
        end

        jBin = local_linearize([jx jy jz], gridShape);
        jHead = spatial.bin_head(jBin);
        if jHead == 0
            continue;
        end

        if jBin == binLin
            i = iHead;
            while i ~= 0
                if activeMask(i)
                    j = spatial.bin_next(i);
                    while j ~= 0
                        if activeMask(j)
                            [keep, dr, r2] = local_pair_geometry(spatial, i, j, cutoff2);
                            if keep
                                nPairs = nPairs + 1;
                                [pairI, pairJ, r2All, drAll] = ...
                                    local_store(nPairs, pairI, pairJ, r2All, drAll, ...
                                                fullToLocal(i), fullToLocal(j), r2, dr, returnDr);
                            end
                        end
                        j = spatial.bin_next(j);
                    end
                end
                i = spatial.bin_next(i);
            end
        else
            i = iHead;
            while i ~= 0
                if activeMask(i)
                    j = jHead;
                    while j ~= 0
                        if activeMask(j)
                            [keep, dr, r2] = local_pair_geometry(spatial, i, j, cutoff2);
                            if keep
                                ii = i;
                                jj = j;
                                li = fullToLocal(ii);
                                lj = fullToLocal(jj);
                                if li < lj
                                    nPairs = nPairs + 1;
                                    [pairI, pairJ, r2All, drAll] = ...
                                        local_store(nPairs, pairI, pairJ, r2All, drAll, ...
                                                    li, lj, r2, dr, returnDr);
                                elseif lj < li
                                    nPairs = nPairs + 1;
                                    [pairI, pairJ, r2All, drAll] = ...
                                        local_store(nPairs, pairI, pairJ, r2All, drAll, ...
                                                    lj, li, r2, -dr, returnDr);
                                end
                            end
                        end
                        j = spatial.bin_next(j);
                    end
                end
                i = spatial.bin_next(i);
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

function [keep, dr, r2] = local_pair_geometry(spatial, i, j, cutoff2)
if spatial.isPeriodic
    fi = spatial.frac(i, :);
    fj = spatial.frac(j, :);
    df = fj - fi;
    df = df - round(df);
    dr = (spatial.cell * df(:)).';
else
    dr = spatial.pos(j, :) - spatial.pos(i, :);
end
r2 = dot(dr, dr);
keep = (r2 <= cutoff2);
end

function [pairI, pairJ, r2All, drAll] = local_store(nPairs, pairI, pairJ, r2All, drAll, i, j, r2, dr, returnDr)
if nPairs > numel(pairI)
    growTo = max(2 * numel(pairI), 1024);
    pairI(growTo, 1) = 0;
    pairJ(growTo, 1) = 0;
    r2All(growTo, 1) = 0;
    if returnDr
        drAll(growTo, 3) = 0;
    end
end
pairI(nPairs) = i;
pairJ(nPairs) = j;
r2All(nPairs) = r2;
if returnDr
    drAll(nPairs, :) = dr;
end
end

function lin = local_linearize(ijk, gridShape)
lin = ijk(1) + (ijk(2)-1)*gridShape(1) + (ijk(3)-1)*gridShape(1)*gridShape(2);
end

function [ix, iy, iz] = local_unlinearize(lin, gridShape)
lin0 = lin - 1;
nxy = gridShape(1) * gridShape(2);
iz = floor(lin0 / nxy) + 1;
rem1 = lin0 - (iz - 1) * nxy;
iy = floor(rem1 / gridShape(1)) + 1;
ix = rem1 - (iy - 1) * gridShape(1) + 1;
end

function value = local_get_opt(s, name, defaultValue)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end