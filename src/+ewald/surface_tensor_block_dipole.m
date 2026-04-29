function Tsurf = surface_tensor_block_dipole(latOrSysOrH, boundary)
%SURFACE_TENSOR_BLOCK_DIPOLE Surface block for dipole FIELD operator.
%
% Tsurf = ewald.surface_tensor_block_dipole(latOrSysOrH)
% Tsurf = ewald.surface_tensor_block_dipole(latOrSysOrH, boundary)
%
% Project lattice convention
% --------------------------
% Polarize uses direct lattice vectors as ROWS:
%
%   r_cart = f_frac * H
%
% This function accepts any geom.get_lattice-compatible input:
%   - raw 3x3 row-lattice H
%   - sys/polsys struct with lattice or super_lattice
%   - lattice struct with fields H, G, and volume
%
% Inputs
%   latOrSysOrH  lattice/system input in project row-lattice convention
%   boundary     'tinfoil' or 'vacuum', default 'tinfoil'
%
% Output
%   Tsurf        3x3 block added between every site pair in the dipole
%                FIELD operator Tpol.
%
% Notes
% -----
% This is the FIELD-operator convention used by the SCF equation
%
%   (I - A*Tpol) mu = A*Eext
%
% For tinfoil boundary conditions:
%
%   Tsurf = 0
%
% For vacuum boundary conditions, the surface contribution to the electric
% field from the total cell dipole M is
%
%   E_surf = -(4*pi/(3V)) * M
%
% Therefore the same pair block added to Tpol is
%
%   Tsurf = -(4*pi/(3V)) * I
%
% This sign is opposite to the most common positive surface ENERGY term.

if nargin < 2 || isempty(boundary)
    boundary = 'tinfoil';
end

if ~(ischar(boundary) || isstring(boundary))
    error('ewald:surface_tensor_block_dipole:BadBoundary', ...
        'boundary must be a character vector or string scalar.');
end

boundary = lower(char(string(boundary)));

lat = local_as_lattice(latOrSysOrH);

if isfield(lat, 'volume') && ~isempty(lat.volume)
    V = lat.volume;
else
    V = abs(det(lat.H));
end

validateattributes(V, {'numeric'}, ...
    {'scalar','real','finite','positive'}, ...
    mfilename, 'cell volume');

if V < 1e-14
    error('ewald:surface_tensor_block_dipole:SingularCell', ...
        'Cell volume is singular or nearly singular.');
end

switch boundary
    case 'tinfoil'
        Tsurf = zeros(3, 3);

    case 'vacuum'
        Tsurf = -(4 * pi / (3 * V)) * eye(3);

    otherwise
        error('ewald:surface_tensor_block_dipole:UnknownBoundary', ...
            'boundary must be ''tinfoil'' or ''vacuum''.');
end
end

function lat = local_as_lattice(x)
if isstruct(x) && isfield(x, 'H') && isfield(x, 'G')
    lat = x;
else
    lat = geom.get_lattice(x);
end

if ~isfield(lat, 'H')
    error('ewald:surface_tensor_block_dipole:BadLattice', ...
        'Input must be compatible with geom.get_lattice and provide H.');
end
end