function test_builder_molecule_indexing()
%TEST_BUILDER_MOLECULE_INDEXING Verify molecule table / site lookup consistency.
%
% This checks the downstream contract used by neighbor selection and charge
% assignment:
%
%   molecule ID -> builder.site_indices_for_molecule -> site set

[filename, tmpDir] = local_write_propene_poscar();
cleanup = onCleanup(@() local_cleanup(tmpDir)); %#ok<NASGU>

crystal = io.import_contcar_as_crystal(filename, ...
    'BondScale', 1.20, ...
    'SortMolecules', false);

model = local_make_model();

reps = [3 2 1];

sys = builder.make_crystal_system(crystal, model, struct( ...
    'supercell_size', reps, ...
    'bondScale', 1.20, ...
    'verbose', false));

nMolExpected = prod(reps) * crystal.nBaseMols;

T = sys.molecule_table;

assert(numel(T.molecule_id) == nMolExpected, ...
    'Unexpected number of reconstructed molecule components.');

assert(isequal(T.molecule_id(:), T.unique_mol_id(:)), ...
    'molecule_id and unique_mol_id aliases should match.');

assert(numel(unique(sys.site_mol_id)) == nMolExpected, ...
    'sys.site_mol_id should contain one label per molecule component.');

% Every molecule component in this non-boundary case should contain 9 sites
% and table indexing should agree with site_indices_for_molecule.
for row = 1:nMolExpected
    molID = T.molecule_id(row);

    idx = builder.site_indices_for_molecule(sys, molID);
    tableIdx = T.site_indices{row};

    assert(isequal(idx(:), tableIdx(:)), ...
        'molecule_table.site_indices should match site_indices_for_molecule.');

    assert(numel(idx) == crystal.nSites, ...
        'Each propene-like molecule component should contain crystal.nSites sites.');

    assert(all(sys.site_mol_id(idx) == molID), ...
        'site_indices_for_molecule returned sites from the wrong molecule.');

    assert(numel(unique(sys.base_mol_id(idx))) == 1, ...
        'A molecule component in this test should have a single base_mol_id.');

    info = builder.molecule_display_status(sys, molID, 'BondScale', 1.20);

    assert(info.is_complete_in_display, ...
        'Each molecule should be complete in this non-boundary-crossing test.');

    assert(info.n_fragments == 1, ...
        'Each molecule should be one displayed fragment.');

    assert(abs(info.largest_fragment_fraction - 1.0) < 1e-12, ...
        'Largest fragment fraction should be 1 for a complete molecule.');

    assert(T.is_complete_in_display(row), ...
        'molecule_table completeness flag should agree with display status.');
end

completeIDs = builder.complete_molecule_ids(sys);

assert(isequal(sort(completeIDs(:)), sort(T.molecule_id(:))), ...
    'All reconstructed molecules should be complete in this test.');

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

fprintf(fid, 'Propene-like molecule indexing test\n');
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