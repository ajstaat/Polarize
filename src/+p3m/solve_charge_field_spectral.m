function [Ex, Ey, Ez, info] = solve_charge_field_spectral(rho, latOrSysOrH, opts)
%SOLVE_CHARGE_FIELD_SPECTRAL Mesh reciprocal field from charge mesh.
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
% The second input may be a lattice struct from geom.get_lattice, a system
% struct, or a 3x3 H matrix.

    if nargin < 3
        opts = struct();
    end

    if ~isfield(opts, 'alpha') || isempty(opts.alpha)
        error('p3m:solve_charge_field_spectral:MissingAlpha', ...
            'opts.alpha is required.');
    end

    alpha = opts.alpha;

    meshSize = size(rho);
    if numel(meshSize) ~= 3
        error('p3m:solve_charge_field_spectral:BadRho', ...
            'rho must be a 3D array.');
    end

    lat = geom.get_lattice(latOrSysOrH);
    H = lat.H;
    G = lat.G; %#ok<NASGU>

    order = local_get_opt(opts, 'assignment_order', 4);
    deconv = local_get_opt(opts, 'deconvolve_assignment', true);
    deconvFloor = local_get_opt(opts, 'deconvolution_floor', 1e-8);

    derivativeMode = lower(char(string(local_get_opt(opts, ...
        'derivative_mode', 'spectral'))));

    fdStencil = lower(char(string(local_get_opt(opts, ...
        'fd_stencil', 'central2'))));

    influenceMode = lower(char(string(local_get_opt(opts, ...
        'influence_mode', 'ewald'))));

    aliasRange = local_get_opt(opts, 'alias_range', 2);

    if ~ismember(derivativeMode, {'spectral','finite_difference'})
        error('p3m:solve_charge_field_spectral:BadDerivativeMode', ...
            'opts.derivative_mode must be ''spectral'' or ''finite_difference''.');
    end

    if strcmp(derivativeMode, 'finite_difference') && ~strcmp(fdStencil, 'central2')
        error('p3m:solve_charge_field_spectral:BadFDStencil', ...
            'Only opts.fd_stencil = ''central2'' is currently supported.');
    end

    if ~ismember(influenceMode, {'ewald','fd_least_squares','optimized'})
        error('p3m:solve_charge_field_spectral:BadInfluenceMode', ...
            ['opts.influence_mode must be ''ewald'', ' ...
             '''fd_least_squares'', or ''optimized''.']);
    end

    validateattributes(aliasRange, {'numeric'}, ...
        {'scalar','real','finite','nonnegative','integer'}, ...
        mfilename, 'opts.alias_range');
    aliasRange = double(aliasRange);

    kg = p3m.make_kgrid(lat, meshSize);

    Qk = fftn(rho);

    D = local_make_derivative_vectors(kg, H, meshSize, derivativeMode, fdStencil);

    influenceInfo = struct();

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
    end

    mask = kg.k2 > 0;

    Ngrid = prod(meshSize);

    switch derivativeMode
        case 'spectral'
            [Ex, Ey, Ez] = local_spectral_field_from_charge_k( ...
                Qk, influence, kg, mask, Ngrid, meshSize);

        case 'finite_difference'
            [Ex, Ey, Ez, phi] = local_finite_difference_field_from_charge_k( ...
                Qk, influence, H, mask, Ngrid, meshSize, fdStencil); %#ok<ASGLU>
    end

    info = struct();
    info.meshSize = meshSize;
    info.alpha = alpha;
    info.assignment_order = order;
    info.deconvolve_assignment = deconv;
    info.deconvolution_floor = deconvFloor;
    info.derivative_mode = derivativeMode;
    info.fd_stencil = fdStencil;
    info.influence_mode = influenceMode;
    info.alias_range = aliasRange;
    info.total_charge_mesh = sum(rho(:));
    info.volume = kg.volume;
    info.nK = nnz(mask);
    info.influenceInfo = influenceInfo;
    info.lattice_convention = 'H_rows_G_columns_HG_2piI';
end

function D = local_make_derivative_vectors(kg, H, meshSize, derivativeMode, fdStencil)
% Return Cartesian derivative vectors D(k) such that grad phi in Fourier
% space is i D(k) phi_k.
%
% With project convention:
%
%   r = s * H
%
% Fractional gradients satisfy:
%
%   grad_s = H * grad_cart
%
% hence:
%
%   grad_cart = H \ grad_s

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

                    A = inv(H);

                    D.x = A(1,1).*D1 + A(1,2).*D2 + A(1,3).*D3;
                    D.y = A(2,1).*D1 + A(2,2).*D2 + A(2,3).*D3;
                    D.z = A(3,1).*D1 + A(3,2).*D2 + A(3,3).*D3;

                otherwise
                    error('p3m:solve_charge_field_spectral:BadFDStencilInternal', ...
                        'Unexpected finite-difference stencil.');
            end

        otherwise
            error('p3m:solve_charge_field_spectral:BadDerivativeModeInternal', ...
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

function [Ex, Ey, Ez] = local_spectral_field_from_charge_k( ...
    Qk, influence, kg, mask, Ngrid, meshSize)

    Fxk = zeros(meshSize);
    Fyk = zeros(meshSize);
    Fzk = zeros(meshSize);

    Fxk(mask) = Ngrid .* (-1i .* kg.kx(mask)) .* influence(mask) .* Qk(mask);
    Fyk(mask) = Ngrid .* (-1i .* kg.ky(mask)) .* influence(mask) .* Qk(mask);
    Fzk(mask) = Ngrid .* (-1i .* kg.kz(mask)) .* influence(mask) .* Qk(mask);

    Ex = real(ifftn(Fxk));
    Ey = real(ifftn(Fyk));
    Ez = real(ifftn(Fzk));
end

function [Ex, Ey, Ez, phi] = local_finite_difference_field_from_charge_k( ...
    Qk, influence, H, mask, Ngrid, meshSize, fdStencil)

    Phik = zeros(meshSize);
    Phik(mask) = Ngrid .* influence(mask) .* Qk(mask);

    phi = real(ifftn(Phik));

    switch fdStencil
        case 'central2'
            [dphi_ds1, dphi_ds2, dphi_ds3] = local_periodic_central2_frac_gradient(phi, meshSize);

        otherwise
            error('p3m:solve_charge_field_spectral:BadFDStencilInternal', ...
                'Unexpected finite-difference stencil.');
    end

    gradS = [dphi_ds1(:).'; dphi_ds2(:).'; dphi_ds3(:).'];
    gradCart = H \ gradS;

    Ecart = -gradCart.';

    Ex = reshape(Ecart(:,1), meshSize);
    Ey = reshape(Ecart(:,2), meshSize);
    Ez = reshape(Ecart(:,3), meshSize);
end

function [df1, df2, df3] = local_periodic_central2_frac_gradient(phi, meshSize)
    M = double(meshSize(:).');

    df1 = 0.5 * M(1) * (circshift(phi, [-1  0  0]) - circshift(phi, [ 1  0  0]));
    df2 = 0.5 * M(2) * (circshift(phi, [ 0 -1  0]) - circshift(phi, [ 0  1  0]));
    df3 = 0.5 * M(3) * (circshift(phi, [ 0  0 -1]) - circshift(phi, [ 0  0  1]));
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

function val = local_get_opt(s, name, defaultVal)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = defaultVal;
    end
end