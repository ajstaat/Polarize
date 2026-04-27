function test_io_import_contcar_as_crystal()
%TEST_IO_IMPORT_CONTCAR_AS_CRYSTAL Verify VASP -> crystal template import.
%
% This test uses a simple propene-like molecule:
%
%   C1H2 = C2H - C3H3
%
% The bond graph does not encode bond order, but the local degrees should be:
%
%   C1 degree 3 -> C_deg3
%   C2 degree 3 -> C_deg3
%   C3 degree 4 -> C_deg4
%
% and hydrogens should be classified by their bonded carbon.

ANG2BOHR = 1.8897259886;

tmpDir = tempname;
mkdir(tmpDir);
cleanup = onCleanup(@() local_cleanup(tmpDir));

filename = fullfile(tmpDir, 'POSCAR_propene');

% Coordinates are Cartesian Angstrom in a large box.
%
% Atom order:
%   1 C1
%   2 C2
%   3 C3
%   4 H on C1
%   5 H on C1
%   6 H on C2
%   7 H on C3
%   8 H on C3
%   9 H on C3
%
% Geometry chosen to make intended bonds unambiguous with a 1.20 covalent
% radius scale:
%   C-C ~1.45 A
%   C-H ~1.09 A
%   nonbonded H/C distances safely larger than cutoff.

xyz = [
    5.00  5.00  5.00   % C1
    6.45  5.00  5.00   % C2
    7.90  5.00  5.00   % C3
    4.55  5.95  5.00   % H on C1
    4.55  4.05  5.00   % H on C1
    6.45  6.09  5.00   % H on C2
    8.35  5.95  5.00   % H on C3
    8.35  4.05  5.00   % H on C3
    7.90  5.00  6.09   % H on C3
];

fid = fopen(filename, 'w');
assert(fid > 0, 'Failed to open temporary POSCAR file.');

fprintf(fid, 'Propene-like import test\n');
fprintf(fid, '1.0\n');
fprintf(fid, '20.0 0.0 0.0\n');
fprintf(fid, '0.0 20.0 0.0\n');
fprintf(fid, '0.0 0.0 20.0\n');
fprintf(fid, 'C H\n');
fprintf(fid, '3 6\n');
fprintf(fid, 'Cartesian\n');

% VASP groups coordinates by species because we declared C H with counts 3 6.
for i = 1:size(xyz, 1)
    fprintf(fid, '%.10f %.10f %.10f\n', xyz(i,1), xyz(i,2), xyz(i,3));
end

fclose(fid);

crystal = io.import_contcar_as_crystal(filename, ...
    'BondScale', 1.20, ...
    'SortMolecules', false);

assert(isstruct(crystal), ...
    'import_contcar_as_crystal should return a struct.');

assert(strcmp(crystal.comment, 'Propene-like import test'), ...
    'crystal.comment was not propagated correctly.');

assert(strcmp(crystal.source_file, filename), ...
    'crystal.source_file was not propagated correctly.');

assert(crystal.nSites == 9, ...
    'Expected 9 sites.');

assert(crystal.nBaseMols == 1, ...
    'Expected one base molecule.');

assert(isequal(size(crystal.lattice), [3 3]), ...
    'crystal.lattice should be 3 x 3.');

assert(isequal(size(crystal.frac_coords), [9 3]), ...
    'crystal.frac_coords should be 9 x 3.');

assert(isequal(size(crystal.cart_coords), [9 3]), ...
    'crystal.cart_coords should be 9 x 3.');

assert(strcmp(crystal.units.length, 'bohr'), ...
    'crystal.units.length should be bohr.');

assert(strcmp(crystal.units.alpha, 'atomic_unit'), ...
    'crystal.units.alpha should be atomic_unit.');

assert(strcmp(crystal.units.charge, 'elementary_charge'), ...
    'crystal.units.charge should be elementary_charge.');

expectedLatticeBohr = 20.0 * ANG2BOHR * eye(3);

assert(norm(crystal.lattice - expectedLatticeBohr, 'fro') < 1e-9, ...
    'crystal.lattice should be converted from Angstrom to bohr.');

assert(norm(crystal.cart_coords - xyz * ANG2BOHR, 'fro') < 1e-9, ...
    'crystal.cart_coords should be converted from Angstrom to bohr.');

assert(norm(crystal.frac_coords - xyz / 20.0, 'fro') < 1e-12, ...
    'crystal.frac_coords should remain dimensionless fractional coordinates.');

assert(all(crystal.base_mol_id == 1), ...
    'All sites should belong to the single propene-like base molecule.');

expectedType = {
    'C'
    'C'
    'C'
    'H'
    'H'
    'H'
    'H'
    'H'
    'H'
};

assert(isequal(crystal.site_type(:), expectedType), ...
    'Unexpected site_type labels.');

expectedClass = {
    'C_deg3'
    'C_deg3'
    'C_deg4'
    'H_on_C_deg3'
    'H_on_C_deg3'
    'H_on_C_deg3'
    'H_on_C_deg4'
    'H_on_C_deg4'
    'H_on_C_deg4'
};

assert(isequal(crystal.site_class(:), expectedClass), ...
    'Unexpected site_class labels.');

expectedLabel = {
    'C1'
    'C2'
    'C3'
    'H1'
    'H2'
    'H3'
    'H4'
    'H5'
    'H6'
};

assert(isequal(crystal.site_label(:), expectedLabel), ...
    'Unexpected molecule-local site labels.');

end

function local_cleanup(tmpDir)
if exist(tmpDir, 'dir')
    rmdir(tmpDir, 's');
end
end