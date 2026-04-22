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
%          .prefer_reach_one    logical, default true for periodic,
%                               false for nonperiodic
%          .directional_pruning logical, default true for both periodic
%                               and nonperiodic
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

preferReachOne     = local_get_opt(opts, 'prefer_reach_one', []);
directionalPruning = local_get_opt(opts, 'directional_pruning', []);
nonperiodicBinScale = local_get_opt(opts, 'nonperiodic_bin_scale', 0.25);

if isempty(preferReachOne)
    preferReachOne = isPeriodic;
end
if isempty(directionalPruning)
    directionalPruning = true;
end

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
    % -------------------------------------------------------------
    % Periodic branch: fractional coordinates, wrapped to [0,1)
    % -------------------------------------------------------------
    frac = (invCell * pos.').';
    frac = frac - floor(frac);

    L = [norm(cellMat(:,1)), norm(cellMat(:,2)), norm(cellMat(:,3))];

    if isempty(binSizeIn)
        baseBinSize = [cutoff cutoff cutoff];

        if preferReachOne
            candidateFactors = [1.00, 1.02, 1.05, 1.07, 1.10, 1.12, 1.15];
            chosenGridShape = [];

            for f = candidateFactors
                trialBinSize = f * baseBinSize;
                trialGridShape = local_periodic_grid_shape_from_bin_size(L, trialBinSize);
                trialBinSizeCart = L .* (1 ./ trialGridShape);
                trialReach = local_periodic_reach_from_bin_size(cutoff, trialBinSizeCart);

                if all(trialReach == [1 1 1])
                    chosenGridShape = trialGridShape;
                    break;
                end
            end

            if isempty(chosenGridShape)
                gridShape = local_periodic_grid_shape_from_bin_size(L, baseBinSize);
            else
                gridShape = chosenGridShape;
            end
        else
            gridShape = local_periodic_grid_shape_from_bin_size(L, baseBinSize);
        end
    else
        if isscalar(binSizeIn)
            binSizeIn = [binSizeIn binSizeIn binSizeIn];
        end
        validateattributes(binSizeIn, {'double'}, ...
            {'vector', 'numel', 3, 'positive', 'finite', 'real'}, ...
            mfilename, 'opts.bin_size');

        gridShape = local_periodic_grid_shape_from_bin_size(L, binSizeIn(:).');
    end

    binSizeFrac = 1 ./ gridShape;
    binSizeCart = L .* binSizeFrac;

    reach = local_periodic_reach_from_bin_size(cutoff, binSizeCart);

    ijk = floor(frac .* gridShape) + 1;
    ijk = min(max(ijk, 1), gridShape);

    if storeFrac
        spatial.frac = frac;
    else
        spatial.frac = [];
    end

    spatial.bin_size = binSizeFrac;
else
    % -------------------------------------------------------------
    % Nonperiodic branch: Cartesian bounding-box binning
    % -------------------------------------------------------------
    xyzMin = min(pos, [], 1);
    xyzMax = max(pos, [], 1);
    extents = max(xyzMax - xyzMin, 0);

    if isempty(binSizeIn)
        binSize = nonperiodicBinScale * [cutoff cutoff cutoff];
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
    spatial.bin_size = binSize;
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

if directionalPruning
    if isPeriodic
        offsets = local_directional_periodic_offsets(cellMat, gridShape, cutoff, reach);
    else
        offsets = local_directional_nonperiodic_offsets(spatial.bin_size, cutoff, reach);
    end
else
    offsets = local_neighbor_offsets(reach);
end

spatial.grid_shape = gridShape;
spatial.n_bins = nBins;
spatial.point_bin = lin;
spatial.bin_head = binHead;
spatial.bin_next = binNext;
spatial.neighbor_offsets = offsets;
spatial.search_reach = reach;

end

function gridShape = local_periodic_grid_shape_from_bin_size(L, binSize)
gridShape = max(ones(1,3), ceil(L ./ binSize));
end

function reach = local_periodic_reach_from_bin_size(cutoff, binSizeCart)
reach = max([1 1 1], ceil(cutoff ./ binSizeCart));
end

function offsets = local_directional_periodic_offsets(cellMat, gridShape, cutoff, reach)
rx = reach(1);
ry = reach(2);
rz = reach(3);

buf = zeros((2*rx+1)*(2*ry+1)*(2*rz+1), 3);
k = 0;

for dx = -rx:rx
    for dy = -ry:ry
        for dz = -rz:rz
            if ~local_keep_half_offset(dx,dy,dz)
                continue;
            end
            if local_keep_offset_directional_periodic(cellMat, gridShape, [dx dy dz], cutoff)
                k = k + 1;
                buf(k,:) = [dx dy dz];
            end
        end
    end
end

offsets = buf(1:k, :);
end

function keep = local_keep_offset_directional_periodic(cellMat, gridShape, dxyz, cutoff)
nx = gridShape(1);
ny = gridShape(2);
nz = gridShape(3);

hx = 1 / nx;
hy = 1 / ny;
hz = 1 / nz;

ds = [dxyz(1)*hx; dxyz(2)*hy; dxyz(3)*hz];
ds = ds - round(ds);

rc = cellMat * ds;
dc = norm(rc);

if dc == 0
    keep = true;
    return;
end

uhat = rc / dc;
q = cellMat.' * uhat;

pad = abs(q(1))*hx + abs(q(2))*hy + abs(q(3))*hz;
keep = (dc <= cutoff + pad);
end

function offsets = local_directional_nonperiodic_offsets(binSize, cutoff, reach)
rx = reach(1);
ry = reach(2);
rz = reach(3);

bx = binSize(1);
by = binSize(2);
bz = binSize(3);

buf = zeros((2*rx+1)*(2*ry+1)*(2*rz+1), 3);
k = 0;

for dx = -rx:rx
    for dy = -ry:ry
        for dz = -rz:rz
            if ~local_keep_half_offset(dx,dy,dz)
                continue;
            end
            if local_keep_offset_directional_nonperiodic([bx by bz], [dx dy dz], cutoff)
                k = k + 1;
                buf(k,:) = [dx dy dz];
            end
        end
    end
end

offsets = buf(1:k, :);
end

function keep = local_keep_offset_directional_nonperiodic(binSize, dxyz, cutoff)
bx = binSize(1);
by = binSize(2);
bz = binSize(3);

rc = [dxyz(1)*bx; dxyz(2)*by; dxyz(3)*bz];
dc = norm(rc);

if dc == 0
    keep = true;
    return;
end

uhat = rc / dc;

% Difference of two within-bin offsets ranges over [-bx,bx], [-by,by], [-bz,bz]
pad = abs(uhat(1))*bx + abs(uhat(2))*by + abs(uhat(3))*bz;

keep = (dc <= cutoff + pad);
end

function tf = local_keep_half_offset(dx,dy,dz)
tf = true;
if dx < 0
    tf = false;
    return;
end
if dx == 0 && dy < 0
    tf = false;
    return;
end
if dx == 0 && dy == 0 && dz < 0
    tf = false;
    return;
end
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