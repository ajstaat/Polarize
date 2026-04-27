function test_builder_make_crystal_system()
%TEST_BUILDER_MAKE_CRYSTAL_SYSTEM Verify crystal template -> supercell system.
%
% This test uses a non-boundary-crossing propene-like molecule. The builder
% should replicate it, reconstruct molecule components, mark all molecules
% complete, and convert input polarizabilities from Angstrom^3 to atomic
% units.

A3_PER_AU_ALPHA = 0.148184711;

[filename, tmpDir] = local_write_propene_poscar();
cleanup = onCleanup(@() local_cleanup(tmpDir)); %#ok<NASGU>

crystal = io.import_contcar_as_crystal(filename, ...
    'BondScale', 1.20, ...
    'SortMolecules', false);

model = local_make_model();

reps = [2 2 1];

sys = builder.make_crystal_system(crystal, model, struct( ...
    'supercell_size', reps, ...
    'bondScale', 1.20, ...
    'verbose', false));

nExpectedSites = crystal.nSites * prod(reps);
nExpectedMols = crystal.nBaseMols * prod(reps);

assert(isstruct(sys), ...
    'builder.make_crystal_system should return a struct.');

io.assert_atomic_units(sys);

assert(strcmp(sys.units.length, 'bohr'), ...
    'sys.units.length should be bohr.');

assert(strcmp(sys.units.alpha, 'atomic_unit'), ...
    'sys.units.alpha should be atomic_unit.');

assert(strcmp(sys.units.charge, 'elementary_charge'), ...
    'sys.units.charge should be elementary_charge.');

assert(isequal(sys.supercell_size, reps), ...
    'sys.supercell_size was not stored correctly.');

assert(sys.n_sites == nExpectedSites, ...
    'Unexpected number of supercell sites.');

assert(size(sys.site_pos, 1) == nExpectedSites && size(sys.site_pos, 2) == 3, ...
    'sys.site_pos should be nSites x 3.');

assert(size(sys.site_frac, 1) == nExpectedSites && size(sys.site_frac, 2) == 3, ...
    'sys.site_frac should be nSites x 3.');

expectedSuperLattice = [
    reps(1) * crystal.lattice(1, :)
    reps(2) * crystal.lattice(2, :)
    reps(3) * crystal.lattice(3, :)
];

assert(norm(sys.super_lattice - expectedSuperLattice, 'fro') < 1e-12, ...
    'Supercell lattice scaling is incorrect.');

assert(numel(sys.base_mol_id) == nExpectedSites, ...
    'sys.base_mol_id should have one entry per site.');

assert(numel(sys.cell_shift) == 3 * nExpectedSites, ...
    'sys.cell_shift should be nSites x 3.');

assert(numel(sys.unit_site_index) == nExpectedSites, ...
    'sys.unit_site_index should have one entry per site.');

assert(numel(sys.site_mol_id) == nExpectedSites, ...
    'sys.site_mol_id should have one entry per site.');

assert(numel(sys.unique_mol_id) == nExpectedSites, ...
    'sys.unique_mol_id should have one entry per site.');

assert(numel(unique(sys.site_mol_id)) == nExpectedMols, ...
    'Unexpected number of reconstructed molecule components.');

assert(numel(sys.site_type) == nExpectedSites, ...
    'sys.site_type should have one entry per site.');

assert(numel(sys.site_class) == nExpectedSites, ...
    'sys.site_class should have one entry per site.');

assert(numel(sys.site_label) == nExpectedSites, ...
    'sys.site_label should have one entry per site.');

assert(numel(sys.site_is_polarizable) == nExpectedSites, ...
    'sys.site_is_polarizable should have one entry per site.');

assert(numel(sys.site_alpha) == nExpectedSites, ...
    'sys.site_alpha should have one entry per site.');

assert(numel(sys.site_charge) == nExpectedSites, ...
    'sys.site_charge should have one entry per site.');

assert(all(sys.site_charge == 0), ...
    'Initial site charges should be zero.');

assert(all(sys.site_is_polarizable), ...
    'All sites in the propene-like test model should be polarizable.');

% Explicit alpha conversion checks. Model values are Angstrom^3 by default;
% sys.site_alpha must be stored in atomic units.
expectedC = 1.750 / A3_PER_AU_ALPHA;
expectedH = 0.696 / A3_PER_AU_ALPHA;

isC = strcmp(sys.site_class, 'C_deg3') | strcmp(sys.site_class, 'C_deg4');
isH = strcmp(sys.site_class, 'H_on_C_deg3') | strcmp(sys.site_class, 'H_on_C_deg4');

assert(all(abs(sys.site_alpha(isC) - expectedC) < 1e-10), ...
    'Carbon polarizabilities should be converted from Angstrom^3 to atomic units.');

assert(all(abs(sys.site_alpha(isH) - expectedH) < 1e-10), ...
    'Hydrogen polarizabilities should be converted from Angstrom^3 to atomic units.');

T = sys.molecule_table;

assert(isstruct(T), ...
    'sys.molecule_table should be a struct.');

assert(isfield(T, 'molecule_id'), ...
    'molecule_table should contain molecule_id.');

assert(isfield(T, 'unique_mol_id'), ...
    'molecule_table should contain unique_mol_id compatibility alias.');

assert(numel(T.molecule_id) == nExpectedMols, ...
    'molecule_table should contain one row per reconstructed molecule.');

assert(numel(T.site_indices) == nExpectedMols, ...
    'molecule_table.site_indices should have one cell per molecule.');

assert(all(T.n_sites == crystal.nSites), ...
    'Each non-boundary propene molecule should have crystal.nSites sites.');

assert(numel(T.is_complete_in_display) == nExpectedMols, ...
    'molecule_table.is_complete_in_display should have one entry per molecule.');

assert(all(T.is_complete_in_display), ...
    'All molecules in this non-boundary-crossing test should be complete.');

completeIDs = builder.complete_molecule_ids(sys);

assert(numel(completeIDs) == nExpectedMols, ...
    'All molecules should be complete in this test.');

end

function model = local_make_model()
model = struct();
model.thole_a = 0.39;
model.alpha_units = 'angstrom^3';

model.polarizable_classes = {
    'C_deg3'
    'C_deg4'
    'H_on_C_deg3'
    'H_on_C_deg4'
};

model.alpha_by_class = struct();
model.alpha_by_class.C_deg3 = 1.750;
model.alpha_by_class.C_deg4 = 1.750;
model.alpha_by_class.H_on_C_deg3 = 0.696;
model.alpha_by_class.H_on_C_deg4 = 0.696;
end

function [filename, tmpDir] = local_write_propene_poscar()
tmpDir = tempname;
mkdir(tmpDir);

filename = fullfile(tmpDir, 'POSCAR_propene');

xyz = [
    5.00  5.00  5.00
    6.45  5.00  5.00
    7.90  5.00  5.00
    4.55  5.95  5.00
    4.55  4.05  5.00
    6.45  6.09  5.00
    8.35  5.95  5.00
    8.35  4.05  5.00
    7.90  5.00  6.09
];

fid = fopen(filename, 'w');
assert(fid > 0, 'Failed to open temporary POSCAR file.');

fprintf(fid, 'Propene-like builder test\n');
fprintf(fid, '1.0\n');
fprintf(fid, '20.0 0.0 0.0\n');
fprintf(fid, '0.0 20.0 0.0\n');
fprintf(fid, '0.0 0.0 20.0\n');
fprintf(fid, 'C H\n');
fprintf(fid, '3 6\n');
fprintf(fid, 'Cartesian\n');

for i = 1:size(xyz, 1)
    fprintf(fid, '%.10f %.10f %.10f\n', xyz(i,1), xyz(i,2), xyz(i,3));
end

fclose(fid);
end

function local_cleanup(tmpDir)
if exist(tmpDir, 'dir')
    rmdir(tmpDir, 's');
end
end