function kCache = build_periodic_kspace_cache(sys, problem, ewaldParams, opts)
%BUILD_PERIODIC_KSPACE_CACHE Build reciprocal-space Ewald cache / plan.
%
% kCache = geom.build_periodic_kspace_cache(sys, problem, ewaldParams)
% kCache = geom.build_periodic_kspace_cache(sys, problem, ewaldParams, opts)
%
% Project lattice convention
% --------------------------
% Polarize uses direct lattice vectors as ROWS:
%
%   r_cart = f_frac * H
%
% geom.get_lattice returns reciprocal vectors as COLUMNS:
%
%   H * G = 2*pi*I
%
% For integer reciprocal index m_col = [h; k; l],
%
%   k_col = G * m_col
%   k_row = k_col.'
%
% This function stores k-vectors as Nk x 3 Cartesian row vectors.
%
% Inputs
%   sys           canonical polarization system with .site_pos and lattice
%   problem       struct with .activeSites
%   ewaldParams   struct with .alpha, .kcut, optional .boundary
%   opts          optional struct:
%                   .kspace_mode             'auto' | 'full' | 'blocked'
%                   .k_block_size            default 2048
%                   .kspace_memory_limit_gb  default 8
%                   .verbose                 default false
%
% Output
%   kCache        reciprocal-space cache for periodic dipole Ewald apply
%
% Notes
% -----
% This cache is for induced-dipole reciprocal-space Ewald, not fixed-charge
% external fields.
%
% The dipole FIELD-operator reciprocal contribution is represented in
% half-k real form. enumerate_kvecs_from_lattice returns one representative
% from each +/- k pair, so two_pref = 2*pref is stored for apply routines.
%
% The prefactor sign here follows the dipole FIELD-operator convention used
% by the periodic induced-dipole routines:
%
%   E_recip = sum_k two_pref(k) * k * [k dot structure_factor]
%
% with pref = -(4*pi/V) exp(-k^2/(4 alpha^2)) / k^2.

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
if strcmp(kspace_mode, 'chunked')
    kspace_mode = 'blocked';
end
if ~ismember(kspace_mode, {'auto', 'full', 'blocked'})
    error('geom:build_periodic_kspace_cache:BadMode', ...
        'opts.kspace_mode must be ''auto'', ''full'', or ''blocked''.');
end

kspace_memory_limit_gb = 8;
if isfield(opts, 'kspace_memory_limit_gb') && ~isempty(opts.kspace_memory_limit_gb)
    kspace_memory_limit_gb = opts.kspace_memory_limit_gb;
end
validateattributes(kspace_memory_limit_gb, {'numeric'}, ...
    {'scalar', 'positive', 'finite', 'real'}, ...
    mfilename, 'opts.kspace_memory_limit_gb');
kspace_memory_limit_gb = double(kspace_memory_limit_gb);

k_block_size = 2048;
if isfield(opts, 'k_block_size') && ~isempty(opts.k_block_size)
    k_block_size = opts.k_block_size;
end
validateattributes(k_block_size, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite', 'real'}, ...
    mfilename, 'opts.k_block_size');
k_block_size = double(k_block_size);

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
    error('geom:build_periodic_kspace_cache:MissingSitePos', ...
        'sys.site_pos is required and may not be empty.');
end

pos = sys.site_pos;
validateattributes(pos, {'numeric'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
    mfilename, 'sys.site_pos');
pos = double(pos);

nSites = size(pos, 1);

if ~isfield(problem, 'activeSites') || isempty(problem.activeSites)
    error('geom:build_periodic_kspace_cache:MissingActiveSites', ...
        'problem.activeSites is required.');
end

activeSites = problem.activeSites(:);
validateattributes(activeSites, {'numeric'}, ...
    {'vector', 'integer', 'positive', 'finite'}, ...
    mfilename, 'problem.activeSites');
activeSites = double(activeSites(:));

if any(activeSites > nSites)
    error('geom:build_periodic_kspace_cache:ActiveSiteOutOfRange', ...
        'problem.activeSites contains indices outside sys.site_pos.');
end

if numel(unique(activeSites)) ~= numel(activeSites)
    error('geom:build_periodic_kspace_cache:DuplicateActiveSites', ...
        'problem.activeSites must not contain duplicate indices.');
end

nPolSites = numel(activeSites);

if ~isfield(ewaldParams, 'alpha') || isempty(ewaldParams.alpha)
    error('geom:build_periodic_kspace_cache:MissingAlpha', ...
        'ewaldParams.alpha is required.');
end
alpha = ewaldParams.alpha;
validateattributes(alpha, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'ewaldParams.alpha');
alpha = double(alpha);

if ~isfield(ewaldParams, 'kcut') || isempty(ewaldParams.kcut)
    error('geom:build_periodic_kspace_cache:MissingKcut', ...
        'ewaldParams.kcut is required.');
end
kcut = ewaldParams.kcut;
validateattributes(kcut, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'ewaldParams.kcut');
kcut = double(kcut);

boundary = 'tinfoil';
if isfield(ewaldParams, 'boundary') && ~isempty(ewaldParams.boundary)
    if ~(ischar(ewaldParams.boundary) || isstring(ewaldParams.boundary))
        error('geom:build_periodic_kspace_cache:BadBoundary', ...
            'ewaldParams.boundary must be a character vector or string scalar.');
    end
    boundary = lower(char(string(ewaldParams.boundary)));
end

if ~ismember(boundary, {'tinfoil', 'vacuum'})
    error('geom:build_periodic_kspace_cache:BadBoundary', ...
        'ewaldParams.boundary must be ''tinfoil'' or ''vacuum''.');
end

lat = geom.get_lattice(sys);
H = lat.H;
G = lat.G;
V = lat.volume;

validateattributes(H, {'numeric'}, {'size', [3 3], 'finite', 'real'}, ...
    mfilename, 'lattice H');
validateattributes(G, {'numeric'}, {'size', [3 3], 'finite', 'real'}, ...
    mfilename, 'lattice G');
validateattributes(V, {'numeric'}, {'scalar', 'finite', 'real', 'positive'}, ...
    mfilename, 'lattice volume');

[kvecs, kmeta] = ewald.enumerate_kvecs_from_lattice(lat, kcut);
nk = size(kvecs, 1);

fullToActive = zeros(nSites, 1);
fullToActive(activeSites) = 1:nPolSites;

pos_pol = pos(activeSites, :);

estimated_full_bytes = double(nPolSites) * double(nk) * 8 * 3;
memory_limit_bytes = kspace_memory_limit_gb * 1024^3;

if strcmp(kspace_mode, 'auto')
    if estimated_full_bytes > memory_limit_bytes
        storage_mode = 'blocked';
    else
        storage_mode = 'full';
    end
else
    storage_mode = kspace_mode;
end

kCache = struct();

kCache.mode = 'periodic_kspace';
kCache.lattice_convention = 'project_row_H_column_G_HG_2piI';
kCache.convention = kmeta.convention;

kCache.nSites = nSites;
kCache.nPolSites = nPolSites;
kCache.activeSites = activeSites;
kCache.full_to_active = fullToActive;

kCache.lattice = lat;
kCache.H = H;
kCache.G = G;
kCache.V = V;

kCache.alpha = alpha;
kCache.kcut = kcut;
kCache.boundary = boundary;

kCache.active_pos = pos_pol;

kCache.k_block_size = k_block_size;
kCache.estimated_full_bytes = estimated_full_bytes;
kCache.estimated_full_gb = estimated_full_bytes / 1024^3;
kCache.memory_limit_bytes = memory_limit_bytes;
kCache.memory_limit_gb = kspace_memory_limit_gb;

kCache.kspace_mode_requested = kspace_mode;
kCache.storage_mode = storage_mode;
kCache.phase_storage = storage_mode;

kCache.kmeta = kmeta;

if nk == 0
    kCache.kvecs = zeros(0, 3);
    kCache.kvecs_T = zeros(3, 0);
    kCache.k2 = zeros(0, 1);
    kCache.knorm = zeros(0, 1);
    kCache.hkl = zeros(0, 3);

    kCache.pref = zeros(0, 1);
    kCache.two_pref = zeros(0, 1);
    kCache.kk6 = zeros(0, 6);

    kCache.cos_phase = zeros(nPolSites, 0);
    kCache.sin_phase = zeros(nPolSites, 0);

    kCache.num_kvec = 0;
    kCache.hkmax = kmeta.hkmax;

    kCache.num_blocks = 0;
    kCache.block_start = zeros(0, 1);
    kCache.block_end = zeros(0, 1);
    kCache.blocks = repmat(local_empty_block(nPolSites), 0, 1);

    if verbose
        fprintf('build_periodic_kspace_cache: no k-vectors within kcut %.6g\n', kcut);
    end

    return;
end

k2 = kmeta.k2(:);
knorm = kmeta.knorm(:);

% Dipole FIELD-operator reciprocal prefactor.
pref = -(4 * pi / V) * exp(-k2 ./ (4 * alpha^2)) ./ k2;
two_pref = 2 * pref;

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
kCache.kvecs_T = kvecs.';
kCache.k2 = k2;
kCache.knorm = knorm;
kCache.hkl = kmeta.hkl;

kCache.pref = pref;
kCache.two_pref = two_pref;
kCache.kk6 = kk6;

kCache.num_kvec = nk;
kCache.hkmax = kmeta.hkmax;

switch storage_mode
    case 'full'
        phase = pos_pol * kvecs.';
        kCache.cos_phase = cos(phase);
        kCache.sin_phase = sin(phase);

        [blockStart, blockEnd, blocks] = local_make_blocks( ...
            kvecs, pref, two_pref, k_block_size, pos_pol, false);

        kCache.num_blocks = numel(blockStart);
        kCache.block_start = blockStart;
        kCache.block_end = blockEnd;
        kCache.blocks = blocks;

        if verbose
            fprintf(['build_periodic_kspace_cache: FULL mode | nPol=%d | nK=%d | ' ...
                     'estimated full storage = %.3f GB\n'], ...
                nPolSites, nk, estimated_full_bytes / 1024^3);
        end

    case 'blocked'
        kCache.cos_phase = [];
        kCache.sin_phase = [];

        [blockStart, blockEnd, blocks] = local_make_blocks( ...
            kvecs, pref, two_pref, k_block_size, pos_pol, true);

        kCache.num_blocks = numel(blockStart);
        kCache.block_start = blockStart;
        kCache.block_end = blockEnd;
        kCache.blocks = blocks;

        if verbose
            phaseBlockBytes = 0;
            for b = 1:kCache.num_blocks
                phaseBlockBytes = phaseBlockBytes + ...
                    8 * numel(kCache.blocks(b).cos_phase) + ...
                    8 * numel(kCache.blocks(b).sin_phase);
            end

            fprintf(['build_periodic_kspace_cache: BLOCKED mode | nPol=%d | nK=%d | ' ...
                     'estimated full storage = %.3f GB | limit = %.3f GB | ' ...
                     'block=%d | nBlocks=%d | stored block phase tables = %.3f GB\n'], ...
                nPolSites, nk, estimated_full_bytes / 1024^3, ...
                memory_limit_bytes / 1024^3, k_block_size, kCache.num_blocks, ...
                phaseBlockBytes / 1024^3);
        end

    otherwise
        error('geom:build_periodic_kspace_cache:InternalModeError', ...
            'Unexpected storage_mode.');
end
end

function [blockStart, blockEnd, blocks] = local_make_blocks( ...
    kvecs, pref, two_pref, k_block_size, pos_pol, store_phase)

nk = size(kvecs, 1);
nPolSites = size(pos_pol, 1);

if nk == 0
    blockStart = zeros(0, 1);
    blockEnd = zeros(0, 1);
    blocks = repmat(local_empty_block(nPolSites), 0, 1);
    return;
end

blockStart = (1:k_block_size:nk).';
nBlocks = numel(blockStart);
blockEnd = zeros(nBlocks, 1);
blocks = repmat(local_empty_block(nPolSites), nBlocks, 1);

for b = 1:nBlocks
    i0 = blockStart(b);
    i1 = min(i0 + k_block_size - 1, nk);
    idx = i0:i1;

    blockEnd(b) = i1;

    blocks(b).idx = idx;
    blocks(b).kvecs = kvecs(idx, :);
    blocks(b).kvecs_T = kvecs(idx, :).';
    blocks(b).pref = pref(idx);
    blocks(b).two_pref = two_pref(idx);
    blocks(b).nk = numel(idx);

    if store_phase
        phase = pos_pol * blocks(b).kvecs_T;
        blocks(b).cos_phase = cos(phase);
        blocks(b).sin_phase = sin(phase);
    else
        blocks(b).cos_phase = zeros(nPolSites, 0);
        blocks(b).sin_phase = zeros(nPolSites, 0);
    end
end
end

function blk = local_empty_block(nPolSites)
if nargin < 1
    nPolSites = 0;
end

blk = struct( ...
    'idx', zeros(1, 0), ...
    'kvecs', zeros(0, 3), ...
    'kvecs_T', zeros(3, 0), ...
    'pref', zeros(0, 1), ...
    'two_pref', zeros(0, 1), ...
    'nk', 0, ...
    'cos_phase', zeros(nPolSites, 0), ...
    'sin_phase', zeros(nPolSites, 0));
end