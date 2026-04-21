function spatial = build_cell_list(pos, opts)
%BUILD_CELL_LIST Construct reusable cell-list spatial index.
%
% spatial = geom.build_cell_list(pos)
% spatial = geom.build_cell_list(pos, opts)
%
% Inputs
%   pos  : N x 3 Cartesian coordinates
%   opts : optional struct with fields
%          .isPeriodic logical, default false
%          .cell       3 x 3 lattice matrix with lattice vectors as COLUMNS
%                      required when isPeriodic = true
%          .cutoff     positive scalar search cutoff used to size bins
%                      required for cell-list mode
%          .bin_size   optional positive scalar or 1x3 vector; default derived
%                      from cutoff
%          .store_frac logical, default true for periodic
%
% Output
%   spatial : struct describing a cell-list index and compatible with
%             geom.query_pairs_cell_list(...)

narginchk(1, 2);
if nargin < 2 || isempty(opts)
    opts = struct();
end

validateattributes(pos, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
    mfilename, 'pos', 1);

isPeriodic = local_get_opt(opts, 'isPeriodic', false);
cellMat    = local_get_opt(opts, 'cell', []);
cutoff     = local_get_opt(opts, 'cutoff', []);
binSizeIn  = local_get_opt(opts, 'bin_size', []);
storeFrac  = local_get_opt(opts, 'store_frac', true);

if isempty(cutoff) || ~(isscalar(cutoff) && isfinite(cutoff) && cutoff > 0)
    error('geom:build_cell_list:BadCutoff', ...
        'opts.cutoff must be a positive finite scalar for cell-list mode.');
end

if isPeriodic
    validateattributes(cellMat, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
        mfilename, 'opts.cell');
    invCell = inv(cellMat);
else
    if ~isempty(cellMat)
        validateattributes(cellMat, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
            mfilename, 'opts.cell');
    end
    invCell = [];
end

n = size(pos, 1);

spatial = struct();
spatial.backend = 'cell_list';
spatial.method = 'cell_list';
spatial.isPeriodic = isPeriodic;
spatial.pos = pos;
spatial.n = n;
spatial.cell = cellMat;
spatial.invCell = invCell;
spatial.cutoff = cutoff;

if n == 0
    spatial.bin_size = [1 1 1];
    spatial.grid_shape = [1 1 1];
    spatial.n_bins = 1;
    spatial.point_bin = zeros(0,1);
    spatial.bin_head = 0;
    spatial.bin_next = zeros(0,1);
    spatial.neighbor_offsets = zeros(0,3);
    spatial.search_reach = [0 0 0];
    if isPeriodic && storeFrac
        spatial.frac = zeros(0,3);
    else
        spatial.frac = [];
    end
    return;
end

if isPeriodic
    % Convert to fractional coordinates and wrap to [0,1).
    frac = (invCell * pos.').';
    frac = frac - floor(frac);

    if isempty(binSizeIn)
        L = [norm(cellMat(:,1)), norm(cellMat(:,2)), norm(cellMat(:,3))];
        gridShape = max(ones(1,3), ceil(L ./ cutoff));
    else
        if isscalar(binSizeIn)
            binSizeIn = [binSizeIn binSizeIn binSizeIn];
        end
        validateattributes(binSizeIn, {'double'}, ...
            {'vector', 'numel', 3, 'positive', 'finite', 'real'}, ...
            mfilename, 'opts.bin_size');
        L = [norm(cellMat(:,1)), norm(cellMat(:,2)), norm(cellMat(:,3))];
        gridShape = max(ones(1,3), ceil(L ./ binSizeIn(:).'));
    end

    binSizeFrac = 1 ./ gridShape;
    binSizeCart = [norm(cellMat(:,1))*binSizeFrac(1), ...
                   norm(cellMat(:,2))*binSizeFrac(2), ...
                   norm(cellMat(:,3))*binSizeFrac(3)];

    reach = max([1 1 1], ceil(cutoff ./ binSizeCart));

    ijk = floor(frac .* gridShape) + 1;
    ijk = min(max(ijk, 1), gridShape);

    if storeFrac
        spatial.frac = frac;
    else
        spatial.frac = [];
    end
else
    xyzMin = min(pos, [], 1);
    xyzMax = max(pos, [], 1);
    extents = max(xyzMax - xyzMin, 0);

    if isempty(binSizeIn)
        binSize = [cutoff cutoff cutoff];
    else
        if isscalar(binSizeIn)
            binSizeIn = [binSizeIn binSizeIn binSizeIn];
        end
        validateattributes(binSizeIn, {'double'}, ...
            {'vector', 'numel', 3, 'positive', 'finite', 'real'}, ...
            mfilename, 'opts.bin_size');
        binSize = binSizeIn(:).';
    end

    gridShape = max(ones(1,3), floor(extents ./ binSize) + 1);
    binSize = extents ./ gridShape;
    binSize(binSize <= 0) = cutoff;

    reach = max([1 1 1], ceil(cutoff ./ binSize));

    rel = pos - xyzMin;
    ijk = floor(rel ./ binSize) + 1;
    ijk = min(max(ijk, 1), gridShape);

    spatial.frac = [];
    spatial.xyz_min = xyzMin;
end

lin = local_linearize_sub(ijk, gridShape);
nBins = prod(gridShape);

binHead = zeros(nBins, 1);
binNext = zeros(n, 1);

for p = n:-1:1
    b = lin(p);
    binNext(p) = binHead(b);
    binHead(b) = p;
end

offsets = local_neighbor_offsets(reach);

if isPeriodic
    spatial.bin_size = binSizeFrac;
else
    spatial.bin_size = binSize;
end
spatial.grid_shape = gridShape;
spatial.n_bins = nBins;
spatial.point_bin = lin;
spatial.bin_head = binHead;
spatial.bin_next = binNext;
spatial.neighbor_offsets = offsets;
spatial.search_reach = reach;

end

function offsets = local_neighbor_offsets(reach)
rx = reach(1);
ry = reach(2);
rz = reach(3);

buf = zeros((2*rx+1)*(2*ry+1)*(2*rz+1), 3);
k = 0;
for dx = -rx:rx
    for dy = -ry:ry
        for dz = -rz:rz
            if dx < 0
                continue;
            end
            if dx == 0 && dy < 0
                continue;
            end
            if dx == 0 && dy == 0 && dz < 0
                continue;
            end
            k = k + 1;
            buf(k,:) = [dx dy dz];
        end
    end
end
offsets = buf(1:k, :);
end

function lin = local_linearize_sub(ijk, gridShape)
lin = ijk(:,1) + ...
      (ijk(:,2)-1) * gridShape(1) + ...
      (ijk(:,3)-1) * gridShape(1) * gridShape(2);
end

function value = local_get_opt(s, name, defaultValue)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end