function frac_coords = cart_to_frac(cart_coords, lattice)
%CART_TO_FRAC Convert Cartesian coordinates to fractional coordinates.
%
% Inputs
%   cart_coords : N x 3 Cartesian coordinates
%   lattice     : 3 x 3 lattice matrix, rows are lattice vectors
%
% Output
%   frac_coords : N x 3 fractional coordinates

    if size(cart_coords, 2) ~= 3
        error('geom:cart_to_frac:BadCoords', ...
            'cart_coords must be N x 3.');
    end

    if ~isequal(size(lattice), [3, 3])
        error('geom:cart_to_frac:BadLattice', ...
            'lattice must be 3 x 3.');
    end

    frac_coords = cart_coords / lattice;
end