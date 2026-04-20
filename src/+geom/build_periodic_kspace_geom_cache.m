function kGeom = build_periodic_kspace_geom_cache(sys, problem, ewaldParams, opts)
%BUILD_PERIODIC_KSPACE_GEOM_CACHE Build shared reciprocal-space geometry cache.
%
% Optimized version:
%   - avoids duplicate source/target work when masks are identical
%   - stores cos/sin tables, not raw phase tables
%
% Notes
%   - Uses the same half-space k enumeration as the current periodic code.
%   - No physics-specific tensor or charge payload is stored here.

    narginchk(3, 4);

    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('geom:build_periodic_kspace_geom_cache:MissingSitePos', ...
            'sys.site_pos is required.');
    end

    pos = sys.site_pos;
    validateattributes(pos, {'double'}, {'2d','ncols',3,'finite','real'}, ...
        mfilename, 'sys.site_pos');

    nSites = size(pos, 1);

    if ~isfield(ewaldParams, 'alpha') || isempty(ewaldParams.alpha)
        error('geom:build_periodic_kspace_geom_cache:MissingAlpha', ...
            'ewaldParams.alpha is required.');
    end
    alpha = ewaldParams.alpha;
    validateattributes(alpha, {'double'}, {'scalar','real','finite','positive'}, ...
        mfilename, 'ewaldParams.alpha');

    if ~isfield(ewaldParams, 'kcut') || isempty(ewaldParams.kcut)
        error('geom:build_periodic_kspace_geom_cache:MissingKcut', ...
            'ewaldParams.kcut is required.');
    end
    kcut = ewaldParams.kcut;
    validateattributes(kcut, {'double'}, {'scalar','real','finite','positive'}, ...
        mfilename, 'ewaldParams.kcut');

    boundary = 'tinfoil';
    if isfield(ewaldParams, 'boundary') && ~isempty(ewaldParams.boundary)
        if ~(ischar(ewaldParams.boundary) || isstring(ewaldParams.boundary))
            error('geom:build_periodic_kspace_geom_cache:BadBoundary', ...
                'ewaldParams.boundary must be a character vector or string scalar.');
        end
        boundary = lower(char(string(ewaldParams.boundary)));
    end

    if isfield(opts, 'target_mask') && ~isempty(opts.target_mask)
        targetMask = logical(opts.target_mask(:));
    else
        if ~isfield(problem, 'polMask') || isempty(problem.polMask)
            error('geom:build_periodic_kspace_geom_cache:MissingPolMask', ...
                'problem.polMask is required when opts.target_mask is not provided.');
        end
        targetMask = logical(problem.polMask(:));
    end

    if isfield(opts, 'source_mask') && ~isempty(opts.source_mask)
        sourceMask = logical(opts.source_mask(:));
    else
        sourceMask = targetMask;
    end

    if numel(targetMask) ~= nSites
        error('geom:build_periodic_kspace_geom_cache:BadTargetMask', ...
            'target_mask must have length nSites.');
    end
    if numel(sourceMask) ~= nSites
        error('geom:build_periodic_kspace_geom_cache:BadSourceMask', ...
            'source_mask must have length nSites.');
    end

    H = local_get_direct_lattice(sys);
    V = abs(det(H));
    if V <= 1e-14
        error('geom:build_periodic_kspace_geom_cache:SingularCell', ...
            'Direct lattice matrix must have nonzero volume.');
    end

    targetSites = find(targetMask);
    sourceSites = find(sourceMask);
    sameMasks = isequal(targetMask, sourceMask);

    [kvecs, meta] = ewald.enumerate_kvecs_triclinic(H, kcut);
    nk = size(kvecs, 1);

    if nk == 0
        kGeom = struct();
        kGeom.mode = 'periodic_kspace_geom';
        kGeom.nSites = nSites;
        kGeom.H = H;
        kGeom.V = V;
        kGeom.alpha = alpha;
        kGeom.kcut = kcut;
        kGeom.boundary = boundary;
        kGeom.kvecs = zeros(0, 3);
        kGeom.k2 = zeros(0, 1);
        kGeom.knorm = zeros(0, 1);
        kGeom.pref_base = zeros(0, 1);
        kGeom.num_kvec = 0;
        kGeom.hkmax = meta.hkmax;
        kGeom.target_mask = targetMask;
        kGeom.source_mask = sourceMask;
        kGeom.targetSites = targetSites;
        kGeom.sourceSites = sourceSites;
        kGeom.target_cos = zeros(numel(targetSites), 0);
        kGeom.target_sin = zeros(numel(targetSites), 0);
        kGeom.source_cos = zeros(numel(sourceSites), 0);
        kGeom.source_sin = zeros(numel(sourceSites), 0);
        kGeom.same_masks = sameMasks;
        return;
    end

    pref_base = (4 * pi / V) * exp(-meta.k2(:) ./ (4 * alpha^2)) ./ meta.k2(:);

    pos_target = pos(targetSites, :);
    phase_target = pos_target * kvecs.';
    target_cos = cos(phase_target);
    target_sin = sin(phase_target);

    if sameMasks
        source_cos = target_cos;
        source_sin = target_sin;
    else
        pos_source = pos(sourceSites, :);
        phase_source = pos_source * kvecs.';
        source_cos = cos(phase_source);
        source_sin = sin(phase_source);
    end

    kGeom = struct();
    kGeom.mode = 'periodic_kspace_geom';
    kGeom.nSites = nSites;
    kGeom.H = H;
    kGeom.V = V;
    kGeom.alpha = alpha;
    kGeom.kcut = kcut;
    kGeom.boundary = boundary;

    kGeom.kvecs = kvecs;
    kGeom.k2 = meta.k2(:);
    kGeom.knorm = meta.knorm(:);
    kGeom.pref_base = pref_base;
    kGeom.num_kvec = nk;
    kGeom.hkmax = meta.hkmax;

    kGeom.target_mask = targetMask;
    kGeom.source_mask = sourceMask;
    kGeom.targetSites = targetSites;
    kGeom.sourceSites = sourceSites;
    kGeom.same_masks = sameMasks;

    kGeom.target_cos = target_cos;
    kGeom.target_sin = target_sin;
    kGeom.source_cos = source_cos;
    kGeom.source_sin = source_sin;
end

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('geom:build_periodic_kspace_geom_cache:MissingLattice', ...
            'sys.super_lattice or sys.lattice is required.');
    end

    validateattributes(H, {'double'}, {'size',[3,3],'finite','real'}, ...
        mfilename, 'lattice');
end