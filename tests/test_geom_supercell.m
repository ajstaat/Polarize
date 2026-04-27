function test_geom_supercell()
%TEST_GEOM_SUPERCELL Verify basic supercell replication and metadata.
%
% Project convention:
%   H rows are direct lattice vectors
%   cart = frac * H

rng(3);

H = [
    10.0   0.0   0.0
     2.0  11.0   0.0
     1.0   3.0  12.0
];

frac = rand(8, 3);
cart = geom.frac_to_cart(frac, H);

unit = struct();
unit.lattice = H;
unit.frac_coords = frac;
unit.cart_coords = cart;

% Use string metadata deliberately. This tests type-safe replication.
unit.site_label = string("X") + string((1:size(frac, 1)).');

% Include common numeric/logical metadata fields.
unit.mol_id = ones(size(frac, 1), 1);
unit.site_is_polarizable = true(size(frac, 1), 1);
unit.site_charge = zeros(size(frac, 1), 1);

reps = [2 3 1];
sc = geom.build_supercell(unit, reps);

nUnit = size(frac, 1);
nExpected = nUnit * prod(reps);

assert(sc.n_unit_sites == nUnit, ...
    'Incorrect number of unit-cell sites reported.');

assert(sc.n_cells == prod(reps), ...
    'Incorrect number of replicated cells reported.');

assert(sc.n_sites == nExpected, ...
    'Incorrect number of supercell sites reported.');

assert(size(sc.frac_coords, 1) == nExpected, ...
    'Incorrect number of supercell fractional coordinates.');

assert(size(sc.cart_coords, 1) == nExpected, ...
    'Incorrect number of supercell Cartesian coordinates.');

assert(isequal(size(sc.lattice), [3, 3]), ...
    'Supercell lattice should be 3 x 3.');

expectedLattice = [
    reps(1) * H(1, :)
    reps(2) * H(2, :)
    reps(3) * H(3, :)
];

assert(norm(sc.lattice - expectedLattice, 'fro') < 1e-12, ...
    'Supercell lattice scaling is incorrect.');

assert(isstring(sc.site_label), ...
    'String site_label metadata should remain string metadata.');

assert(numel(sc.site_label) == nExpected, ...
    'Incorrect number of replicated site labels.');

assert(numel(sc.mol_id) == nExpected, ...
    'Incorrect number of replicated molecule IDs.');

assert(numel(sc.site_is_polarizable) == nExpected, ...
    'Incorrect number of replicated polarizability flags.');

assert(numel(sc.site_charge) == nExpected, ...
    'Incorrect number of replicated charges.');

% Check that first image exactly matches the unit-cell coordinates.
rows1 = 1:nUnit;

assert(norm(sc.cart_coords(rows1, :) - cart, 'fro') < 1e-12, ...
    'First image Cartesian coordinates should match the unit cell.');

assert(norm(sc.frac_coords(rows1, :) - frac ./ reps, 'fro') < 1e-12, ...
    'First image supercell fractional coordinates are incorrect.');

assert(all(sc.cell_shift(rows1, :) == 0, 'all'), ...
    'First image should have zero cell shift.');

assert(all(sc.unit_site_index(rows1) == (1:nUnit).'), ...
    'unit_site_index should identify source unit-cell sites.');

% Check expected number of unique image shifts.
uniqueShifts = unique(sc.cell_shift, 'rows');

assert(size(uniqueShifts, 1) == prod(reps), ...
    'Incorrect number of unique cell shifts.');

end