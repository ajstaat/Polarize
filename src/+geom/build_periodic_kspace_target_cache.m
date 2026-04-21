function tCache = build_periodic_kspace_target_cache(sys, problem, ewaldParams, opts)
%BUILD_PERIODIC_KSPACE_TARGET_CACHE Build reusable target-side reciprocal cache.
%
% Refactored memory-light version:
%   - stores target positions
%   - stores k-space metadata / prefactors
%   - does NOT store full target_cos / target_sin tables
%
% The calling code should evaluate target phases in k-blocks on demand.

    narginchk(3, 4);

    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('geom:build_periodic_kspace_target_cache:MissingSitePos', ...
            'sys.site_pos is required.');
    end

    pos = sys.site_pos;
    validateattributes(pos, {'double'}, {'2d','ncols',3,'finite','real'}, ...
        mfilename, 'sys.site_pos');

    nSites = size(pos, 1);

    if ~isfield(ewaldParams, 'alpha') || isempty(ewaldParams.alpha)
        error('geom:build_periodic_kspace_target_cache:MissingAlpha', ...
            'ewaldParams.alpha is required.');
    end
    alpha = ewaldParams.alpha;

    if ~isfield(ewaldParams, 'kcut') || isempty(ewaldParams.kcut)
        error('geom:build_periodic_kspace_target_cache:MissingKcut', ...
            'ewaldParams.kcut is required.');
    end
    kcut = ewaldParams.kcut;

    boundary = 'tinfoil';
    if isfield(ewaldParams, 'boundary') && ~isempty(ewaldParams.boundary)
        boundary = lower(char(string(ewaldParams.boundary)));
    end

    if isfield(opts, 'target_mask') && ~isempty(opts.target_mask)
        targetMask = logical(opts.target_mask(:));
    else
        if ~isfield(problem, 'polMask') || isempty(problem.polMask)
            error('geom:build_periodic_kspace_target_cache:MissingPolMask', ...
                'problem.polMask is required when opts.target_mask is not provided.');
        end
        targetMask = logical(problem.polMask(:));
    end

    if numel(targetMask) ~= nSites
        error('geom:build_periodic_kspace_target_cache:BadTargetMask', ...
            'target_mask must have length nSites.');
    end

    % Optional chunk size for on-the-fly phase evaluation
    kBlockSize = 2048;
    if isfield(opts, 'k_block_size') && ~isempty(opts.k_block_size)
        kBlockSize = opts.k_block_size;
    end
    validateattributes(kBlockSize, {'double'}, {'scalar','integer','positive'}, ...
        mfilename, 'opts.k_block_size');

    H = local_get_direct_lattice(sys);
    V = abs(det(H));
    if V <= 1e-14
        error('geom:build_periodic_kspace_target_cache:SingularCell', ...
            'Direct lattice matrix must have nonzero volume.');
    end

    targetSites = find(targetMask);
    pos_target = pos(targetSites, :);

    [kvecs, meta] = ewald.enumerate_kvecs_triclinic(H, kcut);
    nk = size(kvecs, 1);

    pref_base = zeros(nk, 1);
    if nk > 0
        pref_base = (4 * pi / V) * exp(-meta.k2(:) ./ (4 * alpha^2)) ./ meta.k2(:);
    end

    tCache = struct();
    tCache.mode = 'periodic_kspace_target';
    tCache.nSites = nSites;
    tCache.H = H;
    tCache.V = V;
    tCache.alpha = alpha;
    tCache.kcut = kcut;
    tCache.boundary = boundary;

    tCache.kvecs = kvecs;
    tCache.k2 = meta.k2(:);
    tCache.knorm = meta.knorm(:);
    tCache.pref_base = pref_base;
    tCache.num_kvec = nk;
    tCache.hkmax = meta.hkmax;

    tCache.target_mask = targetMask;
    tCache.targetSites = targetSites;
    tCache.target_pos = pos_target;

    tCache.k_block_size = kBlockSize;
end

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('geom:build_periodic_kspace_target_cache:MissingLattice', ...
            'sys.super_lattice or sys.lattice is required.');
    end

    validateattributes(H, {'double'}, {'size',[3,3],'finite','real'}, ...
        mfilename, 'lattice');
end