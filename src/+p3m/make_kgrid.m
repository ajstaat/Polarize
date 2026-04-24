function kg = make_kgrid(latOrSysOrH, meshSize)
%MAKE_KGRID Build FFT reciprocal grid using project lattice convention.
%
% Project convention:
%
%   H has direct lattice vectors as rows:
%
%       r_cart = f_frac * H
%
%   G has reciprocal lattice vectors as columns:
%
%       H * G = 2*pi*I
%
%   For integer FFT mode m_col:
%
%       k_col = G * m_col
%       k_row = k_col.'
%
% Inputs
%   latOrSysOrH  lattice struct from geom.get_lattice, system struct, or H
%   meshSize     [M1 M2 M3]
%
% Output fields:
%   kg.m1,kg.m2,kg.m3  integer FFT modes
%   kg.kx,kg.ky,kg.kz  Cartesian reciprocal components
%   kg.k2              |k|^2
%   kg.volume
%   kg.H
%   kg.G

    lat = local_as_lattice(latOrSysOrH);

    M = double(meshSize(:).');
    if numel(M) ~= 3 || any(M < 1) || any(M ~= round(M))
        error('p3m:make_kgrid:BadMeshSize', ...
            'meshSize must be positive integer [M1 M2 M3].');
    end

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
    kg.H = lat.H;
    kg.G = lat.G;
    kg.convention = 'H_rows_G_columns_HG_2piI';
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

function lat = local_as_lattice(x)
    if isstruct(x) && isfield(x, 'H') && isfield(x, 'G')
        lat = x;
    else
        lat = geom.get_lattice(x);
    end
end