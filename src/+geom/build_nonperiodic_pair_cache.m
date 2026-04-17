function geomCache = build_nonperiodic_pair_cache(sys, scfParams, spatial)
%BUILD_NONPERIODIC_PAIR_CACHE Build cached active-site pair geometry.
%
% geomCache = geom.build_nonperiodic_pair_cache(sys, scfParams)
% geomCache = geom.build_nonperiodic_pair_cache(sys, scfParams, spatial)
%
% Inputs
%   sys : struct with fields
%       .site_pos
%       .site_is_polarizable
%
%   scfParams : struct with fields
%       .rcut        scalar cutoff (required)
%
%   spatial : optional spatial index built on sys.site_pos
%
% Output
%   geomCache : struct with cached pair geometry for active polarizable sites
%
% Fields
%   .active_idx
%   .nActive
%   .pair_i
%   .pair_j
%   .dr
%   .r2
%   .r
%   .inv_r
%   .inv_r3
%   .inv_r5
%   .rcut
%   .space
%
% Notes
%   - Indices pair_i and pair_j are in active-space indexing (1:nActive).
%   - This is for the nonperiodic polarization path.
%   - Geometry is cutoff-pruned and cached once for repeated solver use.

    narginchk(2, 3);

    validate_sys(sys);
    validate_scfParams(scfParams);

    if nargin < 3 || isempty(spatial)
        spatial = geom.build_spatial_index(sys.site_pos, struct('isPeriodic', false));
    end

    activeIdx = find(sys.site_is_polarizable);
    nActive = numel(activeIdx);

    queryOpts = struct();
    queryOpts.return_r = true;
    queryOpts.return_dr = true;
    queryOpts.subset_idx = activeIdx;

    pairs = geom.query_pairs_within_cutoff(spatial, scfParams.rcut, queryOpts);

    r = pairs.r;
    r2 = pairs.r2;

    if any(r2 <= 0)
        error('geom:build_nonperiodic_pair_cache:NonPositiveDistance', ...
            'Encountered non-positive pair distance in active-site cache.');
    end

    inv_r  = 1 ./ r;
    inv_r3 = inv_r ./ r2;
    inv_r5 = inv_r3 ./ r2;

    geomCache = struct();
    geomCache.active_idx = activeIdx;
    geomCache.nActive = nActive;

    geomCache.pair_i = pairs.i;
    geomCache.pair_j = pairs.j;

    geomCache.dr = pairs.dr;
    geomCache.r2 = r2;
    geomCache.r = r;
    geomCache.inv_r = inv_r;
    geomCache.inv_r3 = inv_r3;
    geomCache.inv_r5 = inv_r5;

    geomCache.rcut = scfParams.rcut;
    geomCache.space = 'polarizable_only';
end

function validate_sys(sys)
    required = {'site_pos', 'site_is_polarizable'};
    for k = 1:numel(required)
        if ~isfield(sys, required{k})
            error('geom:build_nonperiodic_pair_cache:MissingSysField', ...
                'sys.%s is required.', required{k});
        end
    end

    validateattributes(sys.site_pos, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
        mfilename, 'sys.site_pos');

    validateattributes(sys.site_is_polarizable, {'logical', 'numeric'}, {'vector'}, ...
        mfilename, 'sys.site_is_polarizable');

    if numel(sys.site_is_polarizable) ~= size(sys.site_pos, 1)
        error('geom:build_nonperiodic_pair_cache:SizeMismatch', ...
            'numel(sys.site_is_polarizable) must equal size(sys.site_pos,1).');
    end
end

function validate_scfParams(scfParams)
    if ~isfield(scfParams, 'rcut') || isempty(scfParams.rcut)
        error('geom:build_nonperiodic_pair_cache:MissingRcut', ...
            'scfParams.rcut is required.');
    end
    validateattributes(scfParams.rcut, {'double'}, {'scalar', 'positive', 'finite', 'real'}, ...
        mfilename, 'scfParams.rcut');
end