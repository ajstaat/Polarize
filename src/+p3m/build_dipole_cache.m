function cache = build_dipole_cache(sys, opts)
%BUILD_DIPOLE_CACHE Build cached P3M data for reciprocal dipole fields.
%
% Project lattice convention:
%
%   H rows are direct lattice vectors:
%
%       r_cart = f_frac * H
%
%   G columns are reciprocal lattice vectors:
%
%       H * G = 2*pi*I
%
% This cache stores:
%   - source/target B-spline stencil indices and weights
%   - reciprocal FFT k-grid
%   - dipole reciprocal kernel
%
% Required opts:
%   opts.ewald.alpha
%   opts.mesh_size
%   opts.assignment_order
%
% Optional opts:
%   opts.ewald.kcut
%   opts.use_kcut_mask           default true if opts.ewald.kcut exists
%   opts.target_mask             default all sites
%   opts.source_mask             default polarizable sites if available, else all sites
%   opts.deconvolve_assignment   default true
%   opts.deconvolution_floor     default 1e-8
%   opts.verbose                 default false

    io.assert_atomic_units(sys);

    if nargin < 2
        opts = struct();
    end

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

    lat = geom.get_lattice(sys);

    alpha = opts.ewald.alpha;

    meshSize = double(opts.mesh_size(:).');
    if numel(meshSize) ~= 3 || any(meshSize < 1) || any(meshSize ~= round(meshSize))
        error('p3m:build_dipole_cache:BadMeshSize', ...
            'opts.mesh_size must be positive integer [M1 M2 M3].');
    end

    order = double(opts.assignment_order);
    validateattributes(order, {'numeric'}, ...
        {'scalar','real','finite','positive','integer'}, ...
        mfilename, 'opts.assignment_order');

    nSites = size(sys.site_pos, 1);

    if isfield(opts, 'target_mask') && ~isempty(opts.target_mask)
        targetMask = logical(opts.target_mask(:));
    else
        targetMask = true(nSites, 1);
    end

    if isfield(opts, 'source_mask') && ~isempty(opts.source_mask)
        sourceMask = logical(opts.source_mask(:));
    elseif isfield(sys, 'site_is_polarizable') && ~isempty(sys.site_is_polarizable)
        sourceMask = logical(sys.site_is_polarizable(:));
    else
        sourceMask = true(nSites, 1);
    end

    if numel(targetMask) ~= nSites || numel(sourceMask) ~= nSites
        error('p3m:build_dipole_cache:BadMask', ...
            'target_mask and source_mask must have length N.');
    end

    deconv = local_get_opt(opts, 'deconvolve_assignment', true);
    deconvFloor = local_get_opt(opts, 'deconvolution_floor', 1e-8);

    if isfield(opts, 'use_kcut_mask') && ~isempty(opts.use_kcut_mask)
        useKcutMask = logical(opts.use_kcut_mask);
    else
        useKcutMask = isfield(opts.ewald, 'kcut') && ~isempty(opts.ewald.kcut);
    end

    if useKcutMask
        if ~isfield(opts.ewald, 'kcut') || isempty(opts.ewald.kcut)
            error('p3m:build_dipole_cache:MissingKcut', ...
                'opts.ewald.kcut is required when use_kcut_mask=true.');
        end
        kcut = opts.ewald.kcut;
    else
        kcut = Inf;
    end

    verbose = local_get_opt(opts, 'verbose', false);

    sourceSites = find(sourceMask);
    targetSites = find(targetMask);

    fracSource = sys.site_pos(sourceSites, :) / lat.H;
    fracTarget = sys.site_pos(targetSites, :) / lat.H;

    tBuild = tic;

    sourceStencil = local_build_bspline_stencil(fracSource, meshSize, order);
    targetStencil = local_build_bspline_stencil(fracTarget, meshSize, order);

    kg = local_make_kgrid(lat, meshSize);

    k2 = kg.k2;
    mask = (k2 > 0) & (k2 <= kcut^2);

    kernel = zeros(meshSize);
    kernel(mask) = -(4*pi ./ (lat.volume .* k2(mask))) .* ...
        exp(-k2(mask) ./ (4 * alpha^2));

    if deconv
        W2 = local_assignment_window_squared(kg, meshSize, order);
        W2 = max(W2, deconvFloor);
        kernel(mask) = kernel(mask) ./ W2(mask);
    end

    buildTime = toc(tBuild);

    cache = struct();

    cache.kind = 'p3m_dipole_cache';
    cache.lattice_convention = 'H_rows_G_columns_HG_2piI';

    cache.lat = lat;
    cache.H = lat.H;
    cache.G = lat.G;
    cache.volume = lat.volume;

    cache.nSites = nSites;
    cache.mesh_size = meshSize;
    cache.assignment_order = order;

    cache.alpha = alpha;
    cache.use_kcut_mask = useKcutMask;
    cache.kcut = kcut;
    cache.deconvolve_assignment = deconv;
    cache.deconvolution_floor = deconvFloor;

    cache.target_mask = targetMask;
    cache.source_mask = sourceMask;
    cache.target_sites = targetSites;
    cache.source_sites = sourceSites;

    cache.source_stencil = sourceStencil;
    cache.target_stencil = targetStencil;

    cache.kx = kg.kx;
    cache.ky = kg.ky;
    cache.kz = kg.kz;
    cache.k2 = kg.k2;
    cache.kernel = kernel;
    cache.mask = mask;
    cache.Ngrid = prod(meshSize);
    cache.nK = nnz(mask);
    cache.nK_total_nonzero = nnz(k2 > 0);

    cache.time_build = buildTime;

    if verbose
        fprintf('p3m.build_dipole_cache:\n');
        fprintf('  mesh_size       = [%d %d %d]\n', meshSize);
        fprintf('  order           = %d\n', order);
        fprintf('  alpha           = %.8g\n', alpha);
        fprintf('  use_kcut_mask   = %d\n', useKcutMask);
        fprintf('  kcut            = %.8g\n', kcut);
        fprintf('  nSources        = %d\n', numel(sourceSites));
        fprintf('  nTargets        = %d\n', numel(targetSites));
        fprintf('  nK              = %d\n', cache.nK);
        fprintf('  build time      = %.6f s\n', buildTime);
    end
end

%% ========================================================================
% Local helpers
%% ========================================================================

function st = local_build_bspline_stencil(frac, meshSize, order)
%LOCAL_BUILD_BSPLINE_STENCIL Build cached tensor-product B-spline stencils.
%
% Production convention:
%
%   u = mod(fracPos(a,:),1) .* M
%
% and p3m.bspline_weights_1d expects u in grid coordinates, not fractional
% coordinates.

    n = size(frac, 1);
    nStencil = order^3;

    M = double(meshSize(:).');

    linearIdx = zeros(n, nStencil);
    weight = zeros(n, nStencil);

    for i = 1:n
        u = mod(frac(i, :), 1) .* M;

        [i1, w1] = p3m.bspline_weights_1d(u(1), M(1), order);
        [i2, w2] = p3m.bspline_weights_1d(u(2), M(2), order);
        [i3, w3] = p3m.bspline_weights_1d(u(3), M(3), order);

        c = 0;
        for a = 1:order
            for b = 1:order
                for c3 = 1:order
                    c = c + 1;
                    linearIdx(i, c) = sub2ind(M, i1(a), i2(b), i3(c3));
                    weight(i, c) = w1(a) * w2(b) * w3(c3);
                end
            end
        end
    end

    st = struct();
    st.linear_idx = linearIdx;
    st.weight = weight;
    st.nSites = n;
    st.nStencil = nStencil;
end

function kg = local_make_kgrid(lat, meshSize)

    M = double(meshSize(:).');

    m1v = local_fft_modes(M(1));
    m2v = local_fft_modes(M(2));
    m3v = local_fft_modes(M(3));

    [m1, m2, m3] = ndgrid(m1v, m2v, m3v);

    G = lat.G;

    kx = G(1,1).*m1 + G(1,2).*m2 + G(1,3).*m3;
    ky = G(2,1).*m1 + G(2,2).*m2 + G(2,3).*m3;
    kz = G(3,1).*m1 + G(3,2).*m2 + G(3,3).*m3;

    kg = struct();
    kg.meshSize = M;
    kg.m1 = m1;
    kg.m2 = m2;
    kg.m3 = m3;
    kg.kx = kx;
    kg.ky = ky;
    kg.kz = kz;
    kg.k2 = kx.^2 + ky.^2 + kz.^2;
    kg.volume = lat.volume;
end

function m = local_fft_modes(M)
    M = double(M);

    if mod(M, 2) == 0
        m = [0:(M/2), (-(M/2-1)):-1];
    else
        m = [0:((M-1)/2), (-((M-1)/2)):-1];
    end

    m = double(m);
end

function W2 = local_assignment_window_squared(kg, meshSize, order)
    M = double(meshSize(:).');

    s1 = local_sinc(kg.m1 ./ M(1)).^order;
    s2 = local_sinc(kg.m2 ./ M(2)).^order;
    s3 = local_sinc(kg.m3 ./ M(3)).^order;

    W = s1 .* s2 .* s3;
    W2 = abs(W).^2;
end

function y = local_sinc(x)
    y = ones(size(x));
    mask = abs(x) > 1e-14;
    y(mask) = sin(pi*x(mask)) ./ (pi*x(mask));
end

function val = local_get_opt(s, name, defaultVal)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = defaultVal;
    end
end