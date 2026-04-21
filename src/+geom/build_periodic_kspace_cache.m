function kCache = build_periodic_kspace_cache(sys, problem, ewaldParams, opts)
%BUILD_PERIODIC_KSPACE_CACHE Build reciprocal-space Ewald cache.
%
% kCache = geom.build_periodic_kspace_cache(sys, problem, ewaldParams)
% kCache = geom.build_periodic_kspace_cache(sys, problem, ewaldParams, opts)
%
% Inputs
%   sys         canonical polarization-system struct in atomic units
%   problem     struct from thole.prepare_scf_problem(...)
%   ewaldParams struct with fields:
%       .alpha     Ewald screening parameter
%       .kcut      reciprocal-space cutoff
%       .boundary  optional, default 'tinfoil'
%
%   opts        optional struct with fields:
%       .kspace_mode            'auto' | 'full' | 'chunked' (default 'auto')
%       .kspace_memory_limit_gb positive scalar, default 8
%       .k_block_size           positive integer, default 2048
%       .verbose                logical, default false
%
% Output
%   kCache      struct with reciprocal-space cache data
%
% Fields
%   .mode               'periodic_kspace'
%   .nSites
%   .nPolSites
%   .activeSites
%   .full_to_active
%   .H                  direct lattice matrix (columns are lattice vectors)
%   .V                  cell volume
%   .alpha
%   .kcut
%   .boundary
%   .kvecs              Nk x 3 reciprocal vectors
%   .k2                 Nk x 1 squared norms
%   .knorm              Nk x 1 norms
%   .pref               Nk x 1 reciprocal-space prefactors
%   .kk6                Nk x 6 compact symmetric kk^T entries:
%                       [kxx kyy kzz kxy kxz kyz]
%   .num_kvec
%   .hkmax
%   .storage_mode       'full' | 'chunked'
%   .phase_storage      same as storage_mode
%   .active_pos         nPolSites x 3 active-site positions
%   .k_block_size       block size for chunked mode
%   .estimated_full_bytes
%   .memory_limit_bytes
%
% Full-mode-only fields
%   .cos_phase          nPolSites x Nk
%   .sin_phase          nPolSites x Nk

    narginchk(3, 4);

    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    validateattributes(problem, {'struct'}, {'scalar'}, mfilename, 'problem', 2);
    validateattributes(ewaldParams, {'struct'}, {'scalar'}, mfilename, 'ewaldParams', 3);

    verbose = false;
    if isfield(opts, 'verbose') && ~isempty(opts.verbose)
        verbose = logical(opts.verbose);
    end

    kspace_mode = 'auto';
    if isfield(opts, 'kspace_mode') && ~isempty(opts.kspace_mode)
        kspace_mode = lower(char(string(opts.kspace_mode)));
    end
    if ~ismember(kspace_mode, {'auto','full','chunked'})
        error('geom:build_periodic_kspace_cache:BadMode', ...
            'opts.kspace_mode must be ''auto'', ''full'', or ''chunked''.');
    end

    kspace_memory_limit_gb = 8;
    if isfield(opts, 'kspace_memory_limit_gb') && ~isempty(opts.kspace_memory_limit_gb)
        kspace_memory_limit_gb = opts.kspace_memory_limit_gb;
    end
    validateattributes(kspace_memory_limit_gb, {'double'}, {'scalar', 'positive', 'finite'}, ...
        mfilename, 'opts.kspace_memory_limit_gb');

    k_block_size = 2048;
    if isfield(opts, 'k_block_size') && ~isempty(opts.k_block_size)
        k_block_size = opts.k_block_size;
    end
    validateattributes(k_block_size, {'double'}, {'scalar', 'integer', 'positive'}, ...
        mfilename, 'opts.k_block_size');

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('geom:build_periodic_kspace_cache:MissingSitePos', ...
            'sys.site_pos is required and may not be empty.');
    end
    pos = sys.site_pos;
    validateattributes(pos, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
        mfilename, 'sys.site_pos');

    nSites = size(pos, 1);

    if ~isfield(problem, 'activeSites') || isempty(problem.activeSites)
        error('geom:build_periodic_kspace_cache:MissingActiveSites', ...
            'problem.activeSites is required.');
    end
    activeSites = problem.activeSites(:);
    nPolSites = numel(activeSites);

    if ~isfield(ewaldParams, 'alpha') || isempty(ewaldParams.alpha)
        error('geom:build_periodic_kspace_cache:MissingAlpha', ...
            'ewaldParams.alpha is required.');
    end
    alpha = ewaldParams.alpha;
    validateattributes(alpha, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'ewaldParams.alpha');

    if ~isfield(ewaldParams, 'kcut') || isempty(ewaldParams.kcut)
        error('geom:build_periodic_kspace_cache:MissingKcut', ...
            'ewaldParams.kcut is required.');
    end
    kcut = ewaldParams.kcut;
    validateattributes(kcut, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'ewaldParams.kcut');

    boundary = 'tinfoil';
    if isfield(ewaldParams, 'boundary') && ~isempty(ewaldParams.boundary)
        if ~(ischar(ewaldParams.boundary) || isstring(ewaldParams.boundary))
            error('geom:build_periodic_kspace_cache:BadBoundary', ...
                'ewaldParams.boundary must be a character vector or string scalar.');
        end
        boundary = lower(char(string(ewaldParams.boundary)));
    end

    H = local_get_direct_lattice(sys);
    V = abs(det(H));
    if V <= 1e-14
        error('geom:build_periodic_kspace_cache:SingularCell', ...
            'Direct lattice matrix must have nonzero volume.');
    end

    [kvecs, meta] = ewald.enumerate_kvecs_triclinic(H, kcut);
    nk = size(kvecs, 1);

    fullToActive = zeros(nSites, 1);
    fullToActive(activeSites) = 1:nPolSites;

    pos_pol = pos(activeSites, :);

    % Conservative estimate: phase + cos + sin arrays during full-table construction.
    estimated_full_bytes = double(nPolSites) * double(nk) * 8 * 3;
    memory_limit_bytes = kspace_memory_limit_gb * 1024^3;

    if strcmp(kspace_mode, 'auto')
        if estimated_full_bytes > memory_limit_bytes
            storage_mode = 'chunked';
        else
            storage_mode = 'full';
        end
    else
        storage_mode = kspace_mode;
    end

    kCache = struct();
    kCache.mode = 'periodic_kspace';
    kCache.nSites = nSites;
    kCache.nPolSites = nPolSites;
    kCache.activeSites = activeSites;
    kCache.full_to_active = fullToActive;
    kCache.H = H;
    kCache.V = V;
    kCache.alpha = alpha;
    kCache.kcut = kcut;
    kCache.boundary = boundary;
    kCache.active_pos = pos_pol;
    kCache.k_block_size = k_block_size;
    kCache.estimated_full_bytes = estimated_full_bytes;
    kCache.memory_limit_bytes = memory_limit_bytes;
    kCache.storage_mode = storage_mode;
    kCache.phase_storage = storage_mode;

    if nk == 0
        kCache.kvecs = zeros(0, 3);
        kCache.k2 = zeros(0, 1);
        kCache.knorm = zeros(0, 1);
        kCache.pref = zeros(0, 1);
        kCache.kk6 = zeros(0, 6);
        kCache.cos_phase = zeros(nPolSites, 0);
        kCache.sin_phase = zeros(nPolSites, 0);
        kCache.num_kvec = 0;
        kCache.hkmax = meta.hkmax;
        return;
    end

    k2 = meta.k2(:);
    knorm = meta.knorm(:);

    pref = -(4 * pi / V) * exp(-k2 ./ (4 * alpha^2)) ./ k2;

    kx = kvecs(:, 1);
    ky = kvecs(:, 2);
    kz = kvecs(:, 3);

    kk6 = zeros(nk, 6);
    kk6(:, 1) = kx .* kx;
    kk6(:, 2) = ky .* ky;
    kk6(:, 3) = kz .* kz;
    kk6(:, 4) = kx .* ky;
    kk6(:, 5) = kx .* kz;
    kk6(:, 6) = ky .* kz;

    kCache.kvecs = kvecs;
    kCache.k2 = k2;
    kCache.knorm = knorm;
    kCache.pref = pref;
    kCache.kk6 = kk6;
    kCache.num_kvec = nk;
    kCache.hkmax = meta.hkmax;

    switch storage_mode
        case 'full'
            phase = pos_pol * kvecs.';
            kCache.cos_phase = cos(phase);
            kCache.sin_phase = sin(phase);

            if verbose
                fprintf(['build_periodic_kspace_cache: FULL mode | nPol=%d | nK=%d | ' ...
                         'estimated full storage = %.3f GB\n'], ...
                    nPolSites, nk, estimated_full_bytes / 1024^3);
            end

        case 'chunked'
            kCache.cos_phase = [];
            kCache.sin_phase = [];

            if verbose
                fprintf(['build_periodic_kspace_cache: CHUNKED mode | nPol=%d | nK=%d | ' ...
                         'estimated full storage = %.3f GB exceeds limit %.3f GB | block=%d\n'], ...
                    nPolSites, nk, estimated_full_bytes / 1024^3, ...
                    memory_limit_bytes / 1024^3, k_block_size);
            end

        otherwise
            error('geom:build_periodic_kspace_cache:InternalModeError', ...
                'Unexpected storage_mode.');
    end
end

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('geom:build_periodic_kspace_cache:MissingLattice', ...
            'sys.super_lattice or sys.lattice is required for periodic reciprocal-space cache.');
    end

    validateattributes(H, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
        mfilename, 'lattice');
end