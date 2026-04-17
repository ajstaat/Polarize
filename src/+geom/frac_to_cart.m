function cart_coords = frac_to_cart(frac_coords, lattice)
%FRAC_TO_CART Convert fractional coordinates to Cartesian coordinates.
%
% Inputs
%   frac_coords : N x 3 fractional coordinates
%   lattice     : 3 x 3 lattice matrix, rows are lattice vectors
%
% Output
%   cart_coords : N x 3 Cartesian coordinates

    if size(frac_coords, 2) ~= 3
        error('geom:frac_to_cart:BadCoords', ...
            'frac_coords must be N x 3.');
    end

    if ~isequal(size(lattice), [3, 3])
        error('geom:frac_to_cart:BadLattice', ...
            'lattice must be 3 x 3.');
    end

    cart_coords = frac_coords * lattice;
end