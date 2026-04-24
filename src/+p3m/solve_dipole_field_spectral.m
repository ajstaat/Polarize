function [Ex, Ey, Ez, info] = solve_dipole_field_spectral(Px, Py, Pz, latOrSysOrH, opts)
%SOLVE_DIPOLE_FIELD_SPECTRAL Reciprocal-space mesh field from dipole mesh.
%
% Project lattice convention:
%
%   H has direct lattice vectors as rows:
%
%       r_cart = f_frac * H
%
%   G has reciprocal lattice vectors as columns:
%
%       H * G = 2*pi*I
%
% Computes the periodic Ewald reciprocal dipole field using spectral
% differentiation.

    if nargin < 5
        opts = struct();
    end

    if ~isfield(opts, 'alpha') || isempty(opts.alpha)
        error('p3m:solve_dipole_field_spectral:MissingAlpha', ...
            'opts.alpha is required.');
    end

    alpha = opts.alpha;

    meshSize = size(Px);
    if numel(meshSize) ~= 3
        error('p3m:solve_dipole_field_spectral:BadMesh', ...
            'Px must be a 3D array.');
    end
    if ~isequal(size(Py), meshSize) || ~isequal(size(Pz), meshSize)
        error('p3m:solve_dipole_field_spectral:BadMeshComponents', ...
            'Px, Py, and Pz must have identical sizes.');
    end

    lat = geom.get_lattice(latOrSysOrH);

    order = local_get_opt(opts, 'assignment_order', 4);
    deconv = local_get_opt(opts, 'deconvolve_assignment', true);
    deconvFloor = local_get_opt(opts, 'deconvolution_floor', 1e-8);

    useKcutMask = local_get_opt(opts, 'use_kcut_mask', false);
    if useKcutMask
        if ~isfield(opts, 'kcut') || isempty(opts.kcut)
            error('p3m:solve_dipole_field_spectral:MissingKcut', ...
                'opts.kcut is required when opts.use_kcut_mask=true.');
        end
        kcut = opts.kcut;
    else
        kcut = Inf;
    end

    kg = p3m.make_kgrid(lat, meshSize);

    Pxk = fftn(Px);
    Pyk = fftn(Py);
    Pzk = fftn(Pz);

    k2 = kg.k2;

    mask = (k2 > 0) & (k2 <= kcut^2);

    kernel = zeros(meshSize);
    kernel(mask) = -(4*pi ./ (kg.volume .* k2(mask))) .* ...
        exp(-k2(mask) ./ (4 * alpha^2));

    if deconv
        W2 = local_assignment_window_squared(kg, meshSize, order);
        W2 = max(W2, deconvFloor);
        kernel(mask) = kernel(mask) ./ W2(mask);
    end

    Ngrid = prod(meshSize);

    PdotK = kg.kx .* Pxk + kg.ky .* Pyk + kg.kz .* Pzk;

    Exk = zeros(meshSize);
    Eyk = zeros(meshSize);
    Ezk = zeros(meshSize);

    Exk(mask) = Ngrid .* kg.kx(mask) .* kernel(mask) .* PdotK(mask);
    Eyk(mask) = Ngrid .* kg.ky(mask) .* kernel(mask) .* PdotK(mask);
    Ezk(mask) = Ngrid .* kg.kz(mask) .* kernel(mask) .* PdotK(mask);

    Ex = real(ifftn(Exk));
    Ey = real(ifftn(Eyk));
    Ez = real(ifftn(Ezk));

    info = struct();
    info.meshSize = meshSize;
    info.alpha = alpha;
    info.assignment_order = order;
    info.deconvolve_assignment = deconv;
    info.deconvolution_floor = deconvFloor;
    info.derivative_mode = 'spectral';
    info.influence_mode = 'ewald';
    info.prefactor_convention = 'full_fft_single_k_minus_4pi_over_V';
    info.use_kcut_mask = useKcutMask;
    info.kcut = kcut;
    info.total_Px_mesh = sum(Px(:));
    info.total_Py_mesh = sum(Py(:));
    info.total_Pz_mesh = sum(Pz(:));
    info.total_M_mesh = [info.total_Px_mesh, info.total_Py_mesh, info.total_Pz_mesh];
    info.volume = kg.volume;
    info.nK = nnz(mask);
    info.nK_total_nonzero = nnz(k2 > 0);
    info.lattice_convention = 'H_rows_G_columns_HG_2piI';
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