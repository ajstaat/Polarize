function kCache = build_periodic_kspace_cache(sys, problem, ewaldParams)
%BUILD_PERIODIC_KSPACE_CACHE Build reciprocal-space Ewald cache.
%
% kCache = ewald.build_periodic_kspace_cache(sys, problem, ewaldParams)
%
% Inputs
%   sys         canonical polarization-system struct in atomic units
%   problem     struct from thole.prepare_scf_problem(...)
%   ewaldParams struct with fields:
%       .alpha     Ewald screening parameter
%       .kcut      reciprocal-space cutoff
%       .boundary  optional, default 'tinfoil'
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
%   .cos_phase          nPolSites x Nk
%   .sin_phase          nPolSites x Nk
%   .num_kvec
%   .hkmax
%
% Notes
%   - This cache is active-space focused, since the periodic dipole operator
%     only acts on polarizable sites in the SCF solve.
%   - The reciprocal-space tensor block is built from:
%         2 * pref(k) * cos(k·(r_i-r_j)) * (k k^T)
%     so the natural cached invariants are:
%         * k-space prefactors
%         * compact kk^T entries
%         * per-site phase tables
%   - We intentionally do NOT store k·r_i explicitly, because cos/sin are the
%     repeatedly used quantities for later assembly and matrix-free applies.

    narginchk(3, 3);

    validateattributes(problem, {'struct'}, {'scalar'}, mfilename, 'problem', 2);
    validateattributes(ewaldParams, {'struct'}, {'scalar'}, mfilename, 'ewaldParams', 3);

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('ewald:build_periodic_kspace_cache:MissingSitePos', ...
            'sys.site_pos is required and may not be empty.');
    end
    pos = sys.site_pos;
    validateattributes(pos, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
        mfilename, 'sys.site_pos');

    nSites = size(pos, 1);

    if ~isfield(problem, 'activeSites') || isempty(problem.activeSites)
        error('ewald:build_periodic_kspace_cache:MissingActiveSites', ...
            'problem.activeSites is required.');
    end
    activeSites = problem.activeSites(:);
    nPolSites = numel(activeSites);

    if ~isfield(ewaldParams, 'alpha') || isempty(ewaldParams.alpha)
        error('ewald:build_periodic_kspace_cache:MissingAlpha', ...
            'ewaldParams.alpha is required.');
    end
    alpha = ewaldParams.alpha;
    validateattributes(alpha, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'ewaldParams.alpha');

    if ~isfield(ewaldParams, 'kcut') || isempty(ewaldParams.kcut)
        error('ewald:build_periodic_kspace_cache:MissingKcut', ...
            'ewaldParams.kcut is required.');
    end
    kcut = ewaldParams.kcut;
    validateattributes(kcut, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'ewaldParams.kcut');

    boundary = 'tinfoil';
    if isfield(ewaldParams, 'boundary') && ~isempty(ewaldParams.boundary)
        if ~(ischar(ewaldParams.boundary) || isstring(ewaldParams.boundary))
            error('ewald:build_periodic_kspace_cache:BadBoundary', ...
                'ewaldParams.boundary must be a character vector or string scalar.');
        end
        boundary = lower(char(string(ewaldParams.boundary)));
    end

    H = local_get_direct_lattice(sys);
    V = abs(det(H));
    if V <= 1e-14
        error('ewald:build_periodic_kspace_cache:SingularCell', ...
            'Direct lattice matrix must have nonzero volume.');
    end

    [kvecs, meta] = ewald.enumerate_kvecs_triclinic(H, kcut);
    nk = size(kvecs, 1);

    fullToActive = zeros(nSites, 1);
    fullToActive(activeSites) = 1:nPolSites;

    if nk == 0
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

    % Reciprocal-space prefactor used in the block sum
    %   Tij_recip = sum_k 2*pref(k)*cos(k·(ri-rj))*(k k^T)
    pref = -(4 * pi / V) * exp(-k2 ./ (4 * alpha^2)) ./ k2;

    % Compact symmetric kk^T cache:
    %   [kxx kyy kzz kxy kxz kyz]
    kx = kvecs(:, 1);
    ky = kvecs(:, 2);
    kz = kvecs(:, 3);

    kk6 = zeros(nk, 6);
    kk6(:, 1) = kx .* kx;   % kxx
    kk6(:, 2) = ky .* ky;   % kyy
    kk6(:, 3) = kz .* kz;   % kzz
    kk6(:, 4) = kx .* ky;   % kxy
    kk6(:, 5) = kx .* kz;   % kxz
    kk6(:, 6) = ky .* kz;   % kyz

    % Active-site phase tables.
    pos_pol = pos(activeSites, :);        % nPolSites x 3
    phase = pos_pol * kvecs.';            % nPolSites x Nk

    cos_phase = cos(phase);
    sin_phase = sin(phase);

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

    kCache.kvecs = kvecs;
    kCache.k2 = k2;
    kCache.knorm = knorm;
    kCache.pref = pref;
    kCache.kk6 = kk6;

    kCache.cos_phase = cos_phase;
    kCache.sin_phase = sin_phase;

    kCache.num_kvec = nk;
    kCache.hkmax = meta.hkmax;
end

function H = local_get_direct_lattice(sys)
%LOCAL_GET_DIRECT_LATTICE Extract direct lattice matrix, columns are lattice vectors.

    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('ewald:build_periodic_kspace_cache:MissingLattice', ...
            'sys.super_lattice or sys.lattice is required for periodic reciprocal-space cache.');
    end

    validateattributes(H, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
        mfilename, 'lattice');
end