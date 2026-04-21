function spatial = build_spatial_index(pos, opts)
%BUILD_SPATIAL_INDEX Construct reusable geometric search object.
%
% spatial = geom.build_spatial_index(pos)
% spatial = geom.build_spatial_index(pos, opts)
%
% Inputs
%   pos  : N x 3 Cartesian coordinates
%   opts : optional struct with fields
%          .isPeriodic logical, default false
%          .cell       3 x 3 lattice matrix, required if isPeriodic = true
%                      for periodic geom routines, lattice vectors are columns
%          .method     'auto'|'cell_list'|'bruteforce', default 'auto'
%          .cutoff     positive scalar, required for cell_list, optional for auto
%          .bin_size   optional scalar or 1x3 vector for cell-list bins
%          .store_frac logical, default true
%
% Output
%   spatial : struct suitable for geom.query_pairs_within_cutoff(...)

narginchk(1, 2);
if nargin < 2 || isempty(opts)
    opts = struct();
end

validateattributes(pos, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
    mfilename, 'pos', 1);

isPeriodic = local_get_opt(opts, 'isPeriodic', false);
cellMat    = local_get_opt(opts, 'cell', []);
method     = local_get_opt(opts, 'method', 'auto');
cutoff     = local_get_opt(opts, 'cutoff', []);
binSize    = local_get_opt(opts, 'bin_size', []);
storeFrac  = local_get_opt(opts, 'store_frac', true);

if isPeriodic
    if isempty(cellMat)
        error('geom:build_spatial_index:MissingCell', ...
            'opts.cell is required when opts.isPeriodic is true.');
    end
    validateattributes(cellMat, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
        mfilename, 'opts.cell');
else
    if ~isempty(cellMat)
        validateattributes(cellMat, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
            mfilename, 'opts.cell');
    end
end

if ~ischar(method) && ~isstring(method)
    error('geom:build_spatial_index:BadMethod', ...
        'opts.method must be a character vector or string scalar.');
end
method = char(string(method));

switch method
    case 'auto'
        if isempty(cutoff) || size(pos,1) <= 128
            methodUse = 'bruteforce';
        else
            methodUse = 'cell_list';
        end
    case {'cell_list', 'bruteforce'}
        methodUse = method;
    otherwise
        error('geom:build_spatial_index:UnsupportedMethod', ...
            'Unsupported method "%s".', method);
end

switch methodUse
    case 'bruteforce'
        spatial = struct();
        spatial.backend = 'bruteforce';
        spatial.method = 'bruteforce';
        spatial.pos = pos;
        spatial.n = size(pos, 1);
        spatial.isPeriodic = isPeriodic;
        spatial.cell = cellMat;
        if isPeriodic
            spatial.invCell = inv(cellMat);
        else
            spatial.invCell = [];
        end
        spatial.cutoff = cutoff;
    case 'cell_list'
        spatial = geom.build_cell_list(pos, struct( ...
            'isPeriodic', isPeriodic, ...
            'cell', cellMat, ...
            'cutoff', cutoff, ...
            'bin_size', binSize, ...
            'store_frac', storeFrac));
    otherwise
        error('geom:build_spatial_index:InternalMethodError', ...
            'Unhandled method "%s".', methodUse);
end

end

function value = local_get_opt(s, name, defaultValue)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end