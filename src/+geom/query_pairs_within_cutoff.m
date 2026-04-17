function pairs = query_pairs_within_cutoff(spatial, cutoff, opts)
%QUERY_PAIRS_WITHIN_CUTOFF Find unordered pairs within a distance cutoff.
%
% pairs = geom.query_pairs_within_cutoff(spatial, cutoff, opts)
%
% Inputs
%   spatial : struct from geom.build_spatial_index
%   cutoff  : scalar cutoff distance
%   opts    : optional struct with fields
%       .return_r      logical, default true
%       .return_dr     logical, default true
%       .subset_idx    vector of indices into spatial.pos, default []
%                      If provided, only those points are queried.
%
% Output
%   pairs : struct with fields
%       .i   npairs x 1 indices into subset space if subset_idx given,
%            otherwise indices into spatial.pos
%       .j   npairs x 1
%       .dr  npairs x 3 displacement vectors (j - i), if requested
%       .r2  npairs x 1 squared distances
%       .r   npairs x 1 distances, if requested
%
% Notes
%   - This first draft uses brute-force enumeration with i < j.
%   - For periodic geometry, minimum-image displacement is used.
%   - The subset_idx option is useful for active-site-only queries.

    narginchk(2, 3);

    if nargin < 3 || isempty(opts)
        opts = struct();
    end

    validateattributes(cutoff, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'cutoff', 2);

    return_r  = get_opt(opts, 'return_r', true);
    return_dr = get_opt(opts, 'return_dr', true);
    subsetIdx = get_opt(opts, 'subset_idx', []);

    if isempty(subsetIdx)
        idxMap = (1:spatial.n).';
    else
        validateattributes(subsetIdx, {'numeric'}, {'vector', 'integer', 'positive'}, ...
            mfilename, 'opts.subset_idx');
        subsetIdx = subsetIdx(:);
        if any(subsetIdx > spatial.n)
            error('geom:query_pairs_within_cutoff:SubsetOutOfRange', ...
                'opts.subset_idx contains indices outside spatial.pos.');
        end
        idxMap = subsetIdx;
    end

    pos = spatial.pos(idxMap, :);
    n = size(pos, 1);
    cutoff2 = cutoff * cutoff;

    maxPairs = n * (n - 1) / 2;

    pairI = zeros(maxPairs, 1);
    pairJ = zeros(maxPairs, 1);
    r2All = zeros(maxPairs, 1);

    if return_dr
        drAll = zeros(maxPairs, 3);
    else
        drAll = zeros(0, 3);
    end

    nPairs = 0;

    switch spatial.method
        case 'bruteforce'
            for i = 1:n-1
                ri = pos(i, :);
                for j = i+1:n
                    dr = pos(j, :) - ri;

                    if spatial.isPeriodic
                        dr = minimum_image_displacement(dr, spatial.cell, spatial.invCell);
                    end

                    r2 = dot(dr, dr);
                    if r2 <= cutoff2
                        nPairs = nPairs + 1;
                        pairI(nPairs) = i;
                        pairJ(nPairs) = j;
                        r2All(nPairs) = r2;
                        if return_dr
                            drAll(nPairs, :) = dr;
                        end
                    end
                end
            end

        otherwise
            error('geom:query_pairs_within_cutoff:UnsupportedMethod', ...
                'Unsupported spatial.method "%s".', spatial.method);
    end

    pairI = pairI(1:nPairs);
    pairJ = pairJ(1:nPairs);
    r2All = r2All(1:nPairs);
    if return_dr
        drAll = drAll(1:nPairs, :);
    end

    pairs = struct();
    pairs.i = pairI;
    pairs.j = pairJ;
    pairs.r2 = r2All;

    if return_dr
        pairs.dr = drAll;
    end
    if return_r
        pairs.r = sqrt(r2All);
    end

    pairs.index_space = 'subset';
    pairs.subset_idx = idxMap;
end

function dr = minimum_image_displacement(dr, cellMat, invCell)
%MINIMUM_IMAGE_DISPLACEMENT Apply minimum-image convention in a general cell.

    frac = (invCell * dr(:)).';
    frac = frac - round(frac);
    dr = (cellMat * frac(:)).';
end

function value = get_opt(s, name, defaultValue)
    if isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = defaultValue;
    end
end