function cache = build_dipole_cache(sys, opts)
%BUILD_DIPOLE_CACHE Build cached P3M reciprocal dipole-field operator.
%
% cache = p3m.build_dipole_cache(sys, opts)
%
% Purpose
% -------
% Precompute the geometry, B-spline stencils, FFT k-grid, derivative
% symbols, and reciprocal influence function needed to apply the P3M
% reciprocal dipole field repeatedly:
%
%   source dipoles -> mesh dipoles -> FFT -> reciprocal kernel -> IFFT
%                  -> target field interpolation
%
% This cache is intended for iterative SCF/operator use.
%
% Important
% ---------
% This cache represents the RECIPROCAL P3M dipole field only.
%
% It does not include:
%   - real-space dipole-dipole interactions
%   - local Thole corrections
%   - analytic dipole self term
%   - surface term
%
% The eventual periodic_p3m operator should combine:
%
%   real-space row cache
% + p3m reciprocal dipole cache
% + analytic self block
% + surface block
%
% Project lattice convention
% --------------------------
% Polarize uses direct lattice vectors as ROWS:
%
%   cart = frac * H
%
% geom.get_lattice returns:
%
%   lat.H = H
%   lat.G such that H * G = 2*pi*I
%
% Required opts fields
%   opts.ewald.alpha
%   opts.mesh_size
%   opts.assignment_order
%
% Optional opts fields
%   opts.target_mask              default sys.site_is_polarizable
%   opts.source_mask              default sys.site_is_polarizable
%   opts.derivative_mode          'spectral' or 'finite_difference',
%                                 default 'spectral'
%   opts.fd_stencil               default 'central2'
%   opts.influence_mode           'ewald', 'fd_least_squares',
%                                 or 'optimized', default 'ewald'
%   opts.deconvolve_assignment    default true
%   opts.deconvolution_floor      default 1e-8
%   opts.alias_range              default 2
%   opts.verbose                  default false
%
% Output
%   cache struct consumed by p3m.apply_dipole_cache.

if nargin < 2 || isempty(opts)
    opts = struct();
end

io.assert_atomic_units(sys);

if ~isfield(opts, 'ewald') || isempty(opts.ewald)
    error('p3m:build_dipole_cache:MissingEwald', ...
        'opts.ewald is required.');
end

if ~isfield(opts.ewald, 'alpha') || isempty(opts.ewald.alpha)
    error('p3m:build_dipole_cache:MissingAlpha', ...
        'opts.ewald.alpha is required.');
end

if ~isfield(opts, 'mesh_size') || isempty(opts.mesh_size)
    error('p3m:build_dipole_cache:MissingMeshSize', ...
        'opts.mesh_size is required.');
end

if ~isfield(opts, 'assignment_order') || isempty(opts.assignment_order)
    error('p3m:build_dipole_cache:MissingAssignmentOrder', ...
        'opts.assignment_order is required.');
end

alpha = opts.ewald.alpha;
validateattributes(alpha, {'numeric'}, ...
    {'scalar','real','finite','positive'}, ...
    mfilename, 'opts.ewald.alpha');
alpha = double(alpha);

meshSize = local_validate_mesh_size(opts.mesh_size);
order = local_validate_order(opts.assignment_order);

derivativeMode = lower(char(string(local_get_opt(opts, ...
    'derivative_mode', 'spectral'))));

fdStencil = lower(char(string(local_get_opt(opts, ...
    'fd_stencil', 'central2'))));

influenceMode = lower(char(string(local_get_opt(opts, ...
    'influence_mode', 'ewald'))));

deconv = logical(local_get_opt(opts, 'deconvolve_assignment', true));

deconvFloor = local_get_opt(opts, 'deconvolution_floor', 1e-8);
validateattributes(deconvFloor, {'numeric'}, ...
    {'scalar','real','finite','positive'}, ...
    mfilename, 'opts.deconvolution_floor');
deconvFloor = double(deconvFloor);

aliasRange = local_get_opt(opts, 'alias_range', 2);
validateattributes(aliasRange, {'numeric'}, ...
    {'scalar','real','finite','nonnegative','integer'}, ...
    mfilename, 'opts.alias_range');
aliasRange = double(aliasRange);

verbose = logical(local_get_opt(opts, 'verbose', false));

if ~ismember(derivativeMode, {'spectral','finite_difference'})
    error('p3m:build_dipole_cache:BadDerivativeMode', ...
        'opts.derivative_mode must be ''spectral'' or ''finite_difference''.');
end

if strcmp(derivativeMode, 'finite_difference') && ~strcmp(fdStencil, 'central2')
    error('p3m:build_dipole_cache:BadFDStencil', ...
        'Only opts.fd_stencil = ''central2'' is currently supported.');
end

if ~ismember(influenceMode, {'ewald','fd_least_squares','optimized'})
    error('p3m:build_dipole_cache:BadInfluenceMode', ...
        ['opts.influence_mode must be ''ewald'', ', ...
         '''fd_least_squares'', or ''optimized''.']);
end

pos = sys.site_pos;
validateattributes(pos, {'numeric'}, ...
    {'2d','ncols',3,'real','finite'}, ...
    mfilename, 'sys.site_pos');

nSites = size(pos, 1);

if isfield(sys, 'site_is_polarizable') && ~isempty(sys.site_is_polarizable)
    defaultMask = logical(sys.site_is_polarizable(:));
else
    defaultMask = true(nSites, 1);
end

targetMask = local_get_mask(opts, 'target_mask', defaultMask);
sourceMask = local_get_mask(opts, 'source_mask', defaultMask);

if numel(targetMask) ~= nSites || numel(sourceMask) ~= nSites
    error('p3m:build_dipole_cache:BadMaskSize', ...
        'target_mask and source_mask must have length N.');
end

targetMask = logical(targetMask(:));
sourceMask = logical(sourceMask(:));

targetSites = find(targetMask);
sourceSites = find(sourceMask);

lat = geom.get_lattice(sys);
H = lat.H;

fracAll = pos / H;

fracSource = fracAll(sourceSites, :);
fracTarget = fracAll(targetSites, :);

if verbose
    fprintf('p3m.build_dipole_cache:\n');
    fprintf('  lattice convention       = project_row_H_column_G_HG_2piI\n');
    fprintf('  nSites                   = %d\n', nSites);
    fprintf('  nSources                 = %d\n', numel(sourceSites));
    fprintf('  nTargets                 = %d\n', numel(targetSites));
    fprintf('  mesh_size                = [%d %d %d]\n', meshSize);
    fprintf('  assignment_order         = %d\n', order);
    fprintf('  alpha                    = %.8g\n', alpha);
    fprintf('  derivative_mode          = %s\n', derivativeMode);
    fprintf('  influence_mode           = %s\n', influenceMode);
end

tStart = tic;

sourceStencil = local_build_stencils(fracSource, meshSize, order);
targetStencil = local_build_stencils(fracTarget, meshSize, order);

kg = p3m.make_kgrid(lat, meshSize);
D = local_make_derivative_vectors(kg, H, meshSize, derivativeMode, fdStencil);

switch influenceMode
    case 'ewald'
        [influence, influenceInfo] = local_ewald_influence( ...
            kg, meshSize, alpha, order, deconv, deconvFloor);

    case 'fd_least_squares'
        [influence, influenceInfo] = local_fd_least_squares_influence( ...
            kg, D, meshSize, alpha, order, deconvFloor);

    case 'optimized'
        [influence, influenceInfo] = local_optimized_alias_influence( ...
            kg, D, lat, meshSize, alpha, order, aliasRange, deconvFloor);

    otherwise
        error('p3m:build_dipole_cache:InternalBadInfluenceMode', ...
            'Unexpected influence_mode.');
end

maskK = kg.k2 > 0;
Ngrid = prod(meshSize);

cache = struct();

cache.mode = 'p3m_dipole_reciprocal_cache';
cache.reciprocal_only = true;

cache.nSites = nSites;
cache.nSources = numel(sourceSites);
cache.nTargets = numel(targetSites);

cache.source_mask = sourceMask;
cache.target_mask = targetMask;
cache.source_sites = sourceSites;
cache.target_sites = targetSites;

cache.source_frac = fracSource;
cache.target_frac = fracTarget;

cache.mesh_size = meshSize;
cache.assignment_order = order;
cache.Ngrid = Ngrid;

cache.alpha = alpha;

cache.derivative_mode = derivativeMode;
cache.fd_stencil = fdStencil;
cache.influence_mode = influenceMode;
cache.deconvolve_assignment = deconv;
cache.deconvolution_floor = deconvFloor;
cache.alias_range = aliasRange;

cache.source_stencil = sourceStencil;
cache.target_stencil = targetStencil;

cache.kgrid = kg;
cache.D = D;
cache.influence = influence;
cache.influenceInfo = influenceInfo;
cache.maskK = maskK;

cache.H = lat.H;
cache.G = lat.G;
cache.volume = lat.volume;
cache.lattice = lat;
cache.lattice_convention = 'project_row_H_column_G_HG_2piI';

cache.nK = nnz(maskK);
cache.estimated_mesh_bytes = local_estimated_mesh_bytes(meshSize);
cache.estimated_mesh_gb = cache.estimated_mesh_bytes / 1024^3;

cache.build_time = toc(tStart);

if verbose
    fprintf('  nK                       = %d\n', cache.nK);
    fprintf('  estimated mesh storage   = %.6f GB\n', cache.estimated_mesh_gb);
    fprintf('  build time               = %.6f s\n', cache.build_time);
end
end

% =========================================================================
% Stencils
% =========================================================================

function stencil = local_build_stencils(fracPos, meshSize, order)
n = size(fracPos, 1);

stencil = struct();

stencil.n = n;
stencil.mesh_size = meshSize;
stencil.order = order;

stencil.i1 = zeros(n, order);
stencil.i2 = zeros(n, order);
stencil.i3 = zeros(n, order);

stencil.w1 = zeros(n, order);
stencil.w2 = zeros(n, order);
stencil.w3 = zeros(n, order);

if n == 0
    return;
end

fracWrapped = mod(double(fracPos), 1);

for a = 1:n
    u = fracWrapped(a, :) .* meshSize;

    [i1, w1] = p3m.bspline_weights_1d(u(1), meshSize(1), order);
    [i2, w2] = p3m.bspline_weights_1d(u(2), meshSize(2), order);
    [i3, w3] = p3m.bspline_weights_1d(u(3), meshSize(3), order);

    stencil.i1(a, :) = i1(:).';
    stencil.i2(a, :) = i2(:).';
    stencil.i3(a, :) = i3(:).';

    stencil.w1(a, :) = w1(:).';
    stencil.w2(a, :) = w2(:).';
    stencil.w3(a, :) = w3(:).';
end
end

% =========================================================================
% Derivative / influence helpers
% =========================================================================

function D = local_make_derivative_vectors(kg, H, meshSize, derivativeMode, fdStencil)
switch derivativeMode
    case 'spectral'
        D.x = kg.kx;
        D.y = kg.ky;
        D.z = kg.kz;

    case 'finite_difference'
        switch fdStencil
            case 'central2'
                M = double(meshSize(:).');

                D1 = M(1) .* sin(2*pi .* kg.m1 ./ M(1));
                D2 = M(2) .* sin(2*pi .* kg.m2 ./ M(2));
                D3 = M(3) .* sin(2*pi .* kg.m3 ./ M(3));

                % Project convention: cart = frac * H
                % grad_cart_col = H \ grad_frac_col
                A = H \ eye(3);

                D.x = A(1,1).*D1 + A(1,2).*D2 + A(1,3).*D3;
                D.y = A(2,1).*D1 + A(2,2).*D2 + A(2,3).*D3;
                D.z = A(3,1).*D1 + A(3,2).*D2 + A(3,3).*D3;

            otherwise
                error('p3m:build_dipole_cache:BadFDStencilInternal', ...
                    'Unexpected finite-difference stencil.');
        end

    otherwise
        error('p3m:build_dipole_cache:BadDerivativeModeInternal', ...
            'Unexpected derivative mode.');
end

D.k2 = D.x.^2 + D.y.^2 + D.z.^2;
end

function [influence, info] = local_ewald_influence( ...
    kg, meshSize, alpha, order, deconv, deconvFloor)

influence = zeros(meshSize);

mask = kg.k2 > 0;

influence(mask) = (4*pi ./ (kg.volume * kg.k2(mask))) .* ...
    exp(-kg.k2(mask) ./ (4 * alpha^2));

if deconv
    W2 = local_assignment_window_squared_from_modes( ...
        kg.m1, kg.m2, kg.m3, meshSize, order);

    W2 = max(W2, deconvFloor);

    influence(mask) = influence(mask) ./ W2(mask);
end

info = struct();
info.mode = 'ewald';
info.nAlias = 1;
end

function [influence, info] = local_fd_least_squares_influence( ...
    kg, D, meshSize, alpha, order, deconvFloor)

influence = zeros(meshSize);

mask = (kg.k2 > 0) & (D.k2 > 0);

R = zeros(meshSize);
R(mask) = (4*pi ./ (kg.volume * kg.k2(mask))) .* ...
    exp(-kg.k2(mask) ./ (4 * alpha^2));

DdotK = D.x .* kg.kx + D.y .* kg.ky + D.z .* kg.kz;

W2 = local_assignment_window_squared_from_modes( ...
    kg.m1, kg.m2, kg.m3, meshSize, order);

W2 = max(W2, deconvFloor);

influence(mask) = R(mask) .* DdotK(mask) ./ D.k2(mask) ./ W2(mask);

info = struct();
info.mode = 'fd_least_squares';
info.nAlias = 1;
end

function [influence, info] = local_optimized_alias_influence( ...
    kg, D, lat, meshSize, alpha, order, aliasRange, deconvFloor)

M = double(meshSize(:).');
G = lat.G;

influence = zeros(meshSize);

D2 = D.k2;
validD = D2 > 0;

S0 = zeros(meshSize);
Sx = zeros(meshSize);
Sy = zeros(meshSize);
Sz = zeros(meshSize);

nAlias = 0;

for a1 = -aliasRange:aliasRange
    for a2 = -aliasRange:aliasRange
        for a3 = -aliasRange:aliasRange
            ma1 = kg.m1 + a1 * M(1);
            ma2 = kg.m2 + a2 * M(2);
            ma3 = kg.m3 + a3 * M(3);

            kx = G(1,1).*ma1 + G(1,2).*ma2 + G(1,3).*ma3;
            ky = G(2,1).*ma1 + G(2,2).*ma2 + G(2,3).*ma3;
            kz = G(3,1).*ma1 + G(3,2).*ma2 + G(3,3).*ma3;

            k2a = kx.^2 + ky.^2 + kz.^2;

            aliasMask = k2a > 0;

            U2 = local_assignment_window_squared_from_modes( ...
                ma1, ma2, ma3, meshSize, order);

            R = zeros(meshSize);
            R(aliasMask) = (4*pi ./ (kg.volume .* k2a(aliasMask))) .* ...
                exp(-k2a(aliasMask) ./ (4 * alpha^2));

            S0 = S0 + U2;
            Sx = Sx + U2 .* R .* kx;
            Sy = Sy + U2 .* R .* ky;
            Sz = Sz + U2 .* R .* kz;

            nAlias = nAlias + 1;
        end
    end
end

numerator = D.x .* Sx + D.y .* Sy + D.z .* Sz;
denom = D2 .* max(S0.^2, deconvFloor);

mask = validD & (S0 > 0);

influence(mask) = numerator(mask) ./ denom(mask);
influence(~isfinite(influence)) = 0.0;

info = struct();
info.mode = 'optimized';
info.alias_range = aliasRange;
info.nAlias = nAlias;
info.minS0 = min(S0(:));
info.maxS0 = max(S0(:));
end

function W2 = local_assignment_window_squared_from_modes(m1, m2, m3, meshSize, order)
M = double(meshSize(:).');

s1 = local_sinc(m1 ./ M(1)).^order;
s2 = local_sinc(m2 ./ M(2)).^order;
s3 = local_sinc(m3 ./ M(3)).^order;

W = s1 .* s2 .* s3;
W2 = abs(W).^2;
end

function y = local_sinc(x)
y = ones(size(x));

mask = abs(x) > 1e-14;
y(mask) = sin(pi*x(mask)) ./ (pi*x(mask));
end

% =========================================================================
% Misc helpers
% =========================================================================

function meshSize = local_validate_mesh_size(meshSize)
validateattributes(meshSize, {'numeric'}, ...
    {'vector','numel',3,'integer','positive','finite'}, ...
    mfilename, 'opts.mesh_size');

meshSize = double(meshSize(:).');
end

function order = local_validate_order(order)
validateattributes(order, {'numeric'}, ...
    {'scalar','integer','positive','finite'}, ...
    mfilename, 'opts.assignment_order');

order = double(order);
end

function mask = local_get_mask(opts, name, defaultMask)
if isfield(opts, name) && ~isempty(opts.(name))
    mask = logical(opts.(name)(:));
else
    mask = logical(defaultMask(:));
end
end

function val = local_get_opt(s, name, defaultVal)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = defaultVal;
end
end

function bytes = local_estimated_mesh_bytes(meshSize)
% Approximate working memory for several real and complex mesh arrays.
n = prod(double(meshSize));
bytes = n * 8 * 12;
end