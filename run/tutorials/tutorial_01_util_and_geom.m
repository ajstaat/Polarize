%% tutorial_01_util_and_geom
% This tutorial demonstrates the lowest-level conventions in Polarize:
%
%   1. vectorizing/unvectorizing N x 3 arrays
%   2. row-vector lattice convention
%   3. fractional <-> Cartesian coordinate conversion
%   4. supercell replication
%
% Project lattice convention:
%
%   H rows are direct lattice vectors.
%   cart = frac * H

clear; clc;

fprintf('\n============================================================\n');
fprintf('Tutorial 01: util and core geometry\n');
fprintf('============================================================\n');

%% 1. Stack and unstack N x 3 vectors

mu = randn(5,3);

v = util.stack_xyz(mu);
mu2 = util.unstack_xyz(v);

fprintf('\nVectorization convention:\n');
fprintf('  size(mu) = [%d %d]\n', size(mu,1), size(mu,2));
fprintf('  size(v)  = [%d %d]\n', size(v,1), size(v,2));
fprintf('  roundtrip error = %.3e\n', norm(mu(:) - mu2(:)));

%% 2. Build a lattice

H = [
    10.0   0.0   0.0
     2.0  11.0   0.0
     1.0   3.0  12.0
];

lat = geom.get_lattice(H);

fprintf('\nLattice convention:\n');
fprintf('  H rows are direct lattice vectors\n');
fprintf('  cart = frac * H\n');
fprintf('  volume = %.6f\n', lat.volume);
fprintf('  reciprocal identity error = %.3e\n', lat.identity_error);

%% 3. Convert fractional to Cartesian and back

frac = rand(8,3);

cart = geom.frac_to_cart(frac, H);
frac2 = geom.cart_to_frac(cart, H);

fprintf('\nCoordinate conversion:\n');
fprintf('  frac/cart roundtrip error = %.3e\n', norm(frac - frac2, 'fro'));

%% 4. Replicate a unit cell

unit = struct();
unit.lattice = H;
unit.frac_coords = frac;
unit.cart_coords = cart;
unit.site_label = string("X") + string((1:size(frac,1)).');
unit.mol_id = ones(size(frac,1), 1);
unit.site_is_polarizable = true(size(frac,1), 1);
unit.site_charge = zeros(size(frac,1), 1);

reps = [2 3 1];

sc = geom.build_supercell(unit, reps);

fprintf('\nSupercell replication:\n');
fprintf('  reps             = [%d %d %d]\n', reps);
fprintf('  unit sites       = %d\n', size(unit.frac_coords, 1));
fprintf('  supercell cells  = %d\n', sc.n_cells);
fprintf('  supercell sites  = %d\n', sc.n_sites);

fprintf('\nTutorial 01 completed successfully.\n');