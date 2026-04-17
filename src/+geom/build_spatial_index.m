function spatial = build_spatial_index(pos, opts)
%BUILD_SPATIAL_INDEX Construct reusable geometric search object.
%
% spatial = geom.build_spatial_index(pos, opts)
%
% Inputs
%   pos   : N x 3 Cartesian coordinates
%   opts  : optional struct with fields
%       .isPeriodic   logical, default false
%       .cell         3x3 cell matrix or [], default []
%       .method       'bruteforce' for now, default 'bruteforce'
%
% Output
%   spatial : struct suitable for geom.query_pairs_within_cutoff
%
% Notes
%   - This first draft defines the interface and stores geometry.
%   - The current backend is brute force.
%   - A future cell-list implementation can preserve this interface.

    narginchk(1, 2);

    if nargin < 2 || isempty(opts)
        opts = struct();
    end

    validateattributes(pos, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
        mfilename, 'pos', 1);

    isPeriodic = get_opt(opts, 'isPeriodic', false);
    cellMat    = get_opt(opts, 'cell', []);
    method     = get_opt(opts, 'method', 'bruteforce');

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

    validMethods = {'bruteforce'};
    if ~any(strcmp(method, validMethods))
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

function value = get_opt(s, name, defaultValue)
    if isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = defaultValue;
    end
end