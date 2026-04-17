function spatial = build_spatial_index(pos, opts)
%BUILD_SPATIAL_INDEX Construct reusable geometric search object.
%
% spatial = geom.build_spatial_index(pos)
% spatial = geom.build_spatial_index(pos, opts)
%
% Inputs
%   pos   : N x 3 Cartesian coordinates
%   opts  : optional struct with fields
%       .isPeriodic   logical, default false
%       .cell         3 x 3 lattice matrix, required if isPeriodic = true
%       .method       currently only 'bruteforce', default 'bruteforce'
%
% Output
%   spatial : struct suitable for geom.query_pairs_within_cutoff(...)
%
% Notes
%   - This first version keeps the interface stable and uses a brute-force
%     backend.
%   - A future cell-list backend can preserve the same API.

    narginchk(1, 2);

    if nargin < 2 || isempty(opts)
        opts = struct();
    end

    validateattributes(pos, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
        mfilename, 'pos', 1);

    isPeriodic = local_get_opt(opts, 'isPeriodic', false);
    cellMat    = local_get_opt(opts, 'cell', []);
    method     = local_get_opt(opts, 'method', 'bruteforce');

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

    if ~strcmp(method, 'bruteforce')
        error('geom:build_spatial_index:UnsupportedMethod', ...
            'Unsupported method "%s".', method);
    end

    spatial = struct();
    spatial.pos        = pos;
    spatial.n          = size(pos, 1);
    spatial.isPeriodic = isPeriodic;
    spatial.cell       = cellMat;
    spatial.method     = method;

    if isPeriodic
        spatial.invCell = inv(cellMat);
    else
        spatial.invCell = [];
    end
end

function value = local_get_opt(s, name, defaultValue)
    if isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = defaultValue;
    end
end