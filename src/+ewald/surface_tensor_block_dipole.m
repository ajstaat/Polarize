function Tsurf = surface_tensor_block_dipole(H, boundary)
%SURFACE_TENSOR_BLOCK_DIPOLE Surface-term tensor block for dipole Ewald sum.
%
% Tsurf = ewald.surface_tensor_block_dipole(H)
% Tsurf = ewald.surface_tensor_block_dipole(H, boundary)
%
% Inputs
%   H         3x3 direct lattice matrix, columns are lattice vectors
%   boundary  'tinfoil' or 'vacuum', default 'tinfoil'
%
% Output
%   Tsurf     3x3 block added between every site pair in the dipole operator
%
% Notes
%   For vacuum boundary conditions, the surface energy is
%       Usurf = (2*pi/(3V)) * |sum_i mu_i|^2
%   which corresponds to adding the same block
%       Tsurf = (4*pi/(3V)) * I
%   between every pair (i,j), since
%       U = 0.5 * sum_{ij} mu_i' Tsurf mu_j
%
%   For tinfoil boundary conditions, the surface term is zero.

    if nargin < 2 || isempty(boundary)
        boundary = 'tinfoil';
    end

    validateattributes(H, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
        mfilename, 'H', 1);

    if ~(ischar(boundary) || isstring(boundary))
        error('ewald:surface_tensor_block_dipole:BadBoundary', ...
            'boundary must be a character vector or string scalar.');
    end
    boundary = lower(char(string(boundary)));

    V = abs(det(H));
    if V < 1e-14
        error('ewald:surface_tensor_block_dipole:SingularCell', ...
            'H is singular or nearly singular.');
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