function test_geom_core()
%TEST_GEOM_CORE Verify coordinate and lattice conventions.
%
% Project convention:
%   H rows are direct lattice vectors
%   cart = frac * H

rng(2);

H = [
    10.0   0.0   0.0
     2.0  11.0   0.0
     1.0   3.0  12.0
];

frac = rand(8, 3);

cart = geom.frac_to_cart(frac, H);
frac2 = geom.cart_to_frac(cart, H);

assert(isequal(size(cart), size(frac)), ...
    'frac_to_cart should preserve the coordinate array shape.');

assert(isequal(size(frac2), size(frac)), ...
    'cart_to_frac should preserve the coordinate array shape.');

assert(norm(frac - frac2, 'fro') < 1e-12, ...
    'frac/cart coordinate conversion did not round-trip.');

lat = geom.get_lattice(H);

assert(isstruct(lat), ...
    'geom.get_lattice should return a struct.');

assert(isfield(lat, 'H'), ...
    'geom.get_lattice output missing H.');

assert(isfield(lat, 'G'), ...
    'geom.get_lattice output missing G.');

assert(isfield(lat, 'volume'), ...
    'geom.get_lattice output missing volume.');

assert(isfield(lat, 'identity_error'), ...
    'geom.get_lattice output missing identity_error.');

assert(isequal(size(lat.H), [3, 3]), ...
    'lat.H should be 3 x 3.');

assert(isequal(size(lat.G), [3, 3]), ...
    'lat.G should be 3 x 3.');

assert(norm(lat.H - H, 'fro') < 1e-12, ...
    'lat.H should match the input direct lattice.');

assert(lat.volume > 0, ...
    'Lattice volume should be positive.');

assert(lat.identity_error < 1e-12, ...
    'Direct/reciprocal lattice identity check failed.');

Lmin = geom.shortest_lattice_translation(H);

assert(isfinite(Lmin) && Lmin > 0, ...
    'shortest_lattice_translation should return a positive finite value.');

assert(abs(Lmin - 10.0) < 1e-12, ...
    'Unexpected shortest lattice translation for test lattice.');

% Also check that a cell-parameter construction round-trips to a sensible
% 3 x 3 lattice. Orthorhombic case should be exact.
Hortho = geom.lattice_vectors_from_cell([4, 5, 6, 90, 90, 90]);

expected = [
    4 0 0
    0 5 0
    0 0 6
];

assert(norm(Hortho - expected, 'fro') < 1e-12, ...
    'Orthorhombic lattice construction failed.');

end