function tCache = build_periodic_kspace_target_cache(sys, problem, ewaldParams, opts)
%BUILD_PERIODIC_KSPACE_TARGET_CACHE Build reusable target-side reciprocal cache.
%
% tCache = geom.build_periodic_kspace_target_cache(sys, problem, ewaldParams)
% tCache = geom.build_periodic_kspace_target_cache(sys, problem, ewaldParams, opts)
%
% Inputs
%   sys         canonical polarization-system struct
%   problem     struct from thole.prepare_scf_problem(...)
%   ewaldParams struct with fields:
%       .alpha
%       .kcut
%       .boundary   optional, default 'tinfoil'
%
%   opts        optional struct with fields:
%       .target_mask   N x 1 logical, optional
%
% Output
%   tCache struct with fields:
%       .mode
%       .nSites
%       .H
%       .V
%       .alpha
%       .kcut
%       .boundary
%       .kvecs
%       .k2
%       .knorm
%       .pref_base
%       .num_kvec
%       .hkmax
%       .target_mask
%       .targetSites
%       .target_cos
%       .target_sin

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

    H = local_get_direct_lattice(sys);
    V = abs(det(H));
    if V <= 1e-14
        error('geom:build_periodic_kspace_target_cache:SingularCell', ...
            'Direct lattice matrix must have nonzero volume.');
    end

    targetSites = find(targetMask);

    [kvecs, meta] = ewald.enumerate_kvecs_triclinic(H, kcut);
    nk = size(kvecs, 1);

    pref_base = zeros(nk, 1);
    target_cos = zeros(numel(targetSites), nk);
    target_sin = zeros(numel(targetSites), nk);

    if nk > 0
        pref_base = (4 * pi / V) * exp(-meta.k2(:) ./ (4 * alpha^2)) ./ meta.k2(:);
        pos_target = pos(targetSites, :);
        phase_target = pos_target * kvecs.';
        target_cos = cos(phase_target);
        target_sin = sin(phase_target);
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
    tCache.target_cos = target_cos;
    tCache.target_sin = target_sin;
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