function kg = make_kgrid(latOrSysOrH, meshSize)
%MAKE_KGRID Build FFT reciprocal grid using the project lattice convention.
%
% kg = p3m.make_kgrid(latOrSysOrH, meshSize)
%
% Project lattice convention
% --------------------------
% Polarize uses direct lattice vectors as ROWS:
%
%   cart = frac * H
%
% geom.get_lattice returns reciprocal vectors as COLUMNS:
%
%   H * G = 2*pi*I
%
% For integer FFT mode m_col = [m1; m2; m3]:
%
%   k_col = G * m_col
%   k_row = k_col.'
%
% Inputs
%   latOrSysOrH   one of:
%                   - lattice struct from geom.get_lattice
%                   - system struct accepted by geom.get_lattice
%                   - 3 x 3 row-lattice matrix H
%
%   meshSize      [M1 M2 M3] positive integer FFT mesh dimensions
%
% Output
%   kg struct with fields:
%     .meshSize
%     .m1,.m2,.m3       integer FFT mode grids
%     .kx,.ky,.kz       Cartesian reciprocal components
%     .k2               |k|^2
%     .volume
%     .H
%     .G
%     .lattice
%     .convention

lat = local_as_lattice(latOrSysOrH);

M = local_validate_mesh_size(meshSize);

G = lat.G;

m1v = local_fft_modes(M(1));
m2v = local_fft_modes(M(2));
m3v = local_fft_modes(M(3));

[m1, m2, m3] = ndgrid(m1v, m2v, m3v);

% k_col = G * [m1;m2;m3]
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
kg.lattice = lat;

kg.convention = 'project_row_H_column_G_HG_2piI';

if isfield(lat, 'convention')
    kg.lattice_convention = lat.convention;
else
    kg.lattice_convention = kg.convention;
end
end

function m = local_fft_modes(M)
M = double(M);

if mod(M, 2) == 0
    % MATLAB FFT ordering for even M:
    %   0, 1, ..., M/2, -(M/2-1), ..., -1
    %
    % The Nyquist mode M/2 is represented once.
    m = [0:(M/2), (-(M/2-1)):-1];
else
    m = [0:((M-1)/2), (-((M-1)/2)):-1];
end

m = double(m);
end

function lat = local_as_lattice(x)
if isstruct(x) && isfield(x, 'H') && isfield(x, 'G') && isfield(x, 'volume')
    lat = x;
else
    lat = geom.get_lattice(x);
end

validateattributes(lat.H, {'numeric'}, ...
    {'size', [3 3], 'real', 'finite'}, ...
    mfilename, 'lat.H');

validateattributes(lat.G, {'numeric'}, ...
    {'size', [3 3], 'real', 'finite'}, ...
    mfilename, 'lat.G');

validateattributes(lat.volume, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'lat.volume');

err = norm(lat.H * lat.G - 2*pi*eye(3), 'fro');
if err > 1e-9
    error('p3m:make_kgrid:BadLatticeConvention', ...
        'Expected H*G = 2*pi*I. ||H*G-2*piI||_F = %.3e.', err);
end
end

function M = local_validate_mesh_size(meshSize)
validateattributes(meshSize, {'numeric'}, ...
    {'vector', 'numel', 3, 'integer', 'positive', 'finite'}, ...
    mfilename, 'meshSize');

M = double(meshSize(:).');
end