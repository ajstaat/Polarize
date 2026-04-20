function Usurf = surface_term_dipole(mu, H, boundary)
%SURFACE_TERM_DIPOLE Surface term for dipole Ewald sum.
%
% Usurf = ewald.surface_term_dipole(mu, H)
% Usurf = ewald.surface_term_dipole(mu, H, boundary)
%
% Inputs
%   mu        N x 3 dipole array
%   H         3x3 direct lattice matrix, columns are lattice vectors
%   boundary  'tinfoil' or 'vacuum', default 'tinfoil'
%
% Output
%   Usurf     scalar surface-energy contribution
%
% Notes
%   - For conducting ('tinfoil') boundary conditions, the surface term is zero.
%   - For vacuum boundary conditions:
%       Usurf = (2*pi / (3*V)) * |sum_i mu_i|^2

    if nargin < 3 || isempty(boundary)
        boundary = 'tinfoil';
    end

    validateattributes(mu, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
        mfilename, 'mu', 1);
    validateattributes(H, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
        mfilename, 'H', 2);

    if ~(ischar(boundary) || isstring(boundary))
        error('ewald:surface_term_dipole:BadBoundary', ...
            'boundary must be a character vector or string scalar.');
    end
    boundary = lower(char(string(boundary)));

    V = abs(det(H));
    if V < 1e-14
        error('ewald:surface_term_dipole:SingularCell', ...
            'H is singular or nearly singular.');
    end

    M = sum(mu, 1);

    switch boundary
        case 'tinfoil'
            Usurf = 0.0;
        case 'vacuum'
            Usurf = (2 * pi / (3 * V)) * dot(M, M);
        otherwise
            error('ewald:surface_term_dipole:UnknownBoundary', ...
                'boundary must be ''tinfoil'' or ''vacuum''.');
    end
end