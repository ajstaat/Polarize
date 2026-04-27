function test_builder_boundary_molecule_reconstruction()
%TEST_BUILDER_BOUNDARY_MOLECULE_RECONSTRUCTION Verify supercell graph reconstruction.
%
% This test constructs a propene-like molecule split across the unit-cell x
% boundary. IO should identify it as one base molecule using the unit-cell
% PBC graph.
%
% Builder then replicates the unit cell and uses the full-supercell PBC
% graph to reconstruct molecular components. In a sufficiently replicated
% supercell, complete displayed molecules should exist even though the
% original unit-cell template is wrapped.

[filename, tmpDir] = local_write_boundary_propene_poscar();
cleanup = onCleanup(@() local_cleanup(tmpDir)); %#ok<NASGU>

crystal = io.import_contcar_as_crystal(filename, ...
    'BondScale', 1.20, ...
    'SortMolecules', false);

assert(crystal.nSites == 9, ...
    'Expected 9 sites in the boundary-crossing propene-like molecule.');

assert(crystal.nBaseMols == 1, ...
    'IO should identify the boundary-crossing molecule as one base molecule.');

assert(all(crystal.base_mol_id == 1), ...
    'All sites should belong to one base molecule after unit-cell PBC import.');

model = local_make_model();

% A [4 1 1] supercell should contain enough replicated images for the
% supercell PBC graph + display completeness checker to find complete
% interior molecule components.
reps = [4 1 1];

sys = builder.make_crystal_system(crystal, model, struct( ...
    'supercell_size', reps, ...
    'bondScale', 1.20, ...
    'verbose', false));

assert(sys.n_sites == crystal.nSites * prod(reps), ...
    'Unexpected number of supercell sites.');

T = sys.molecule_table;

assert(isfield(T, 'molecule_id'), ...
    'molecule_table should contain molecule_id.');

assert(isfield(T, 'is_complete_in_display'), ...
    'molecule_table should contain is_complete_in_display.');

% The exact number of connected components can include boundary-wrapping
% components/fragments, but complete components should exist.
completeIDs = builder.complete_molecule_ids(sys);

nExpectedComponents = prod(reps) * crystal.nBaseMols;
nExpectedComplete = reps(1) - 1;  % for this 1D x-boundary-crossing setup

assert(numel(T.molecule_id) == nExpectedComponents, ...
    'Expected one reconstructed PBC molecule component per replicated unit cell.');

assert(numel(completeIDs) == nExpectedComplete, ...
    'Expected exactly three complete interior molecules in the [4 1 1] boundary-crossing test.');

assert(nnz(~T.is_complete_in_display) == 1, ...
    'Expected exactly one incomplete boundary-wrapping molecule in the [4 1 1] test.');

% Every complete molecule should have the full propene site count and be one
% displayed fragment.
for k = 1:numel(completeIDs)
    molID = completeIDs(k);
    idx = builder.site_indices_for_molecule(sys, molID);

    assert(numel(idx) == crystal.nSites, ...
        'Each complete propene molecule should have the full 9-site count.');

    info = builder.molecule_display_status(sys, molID, 'BondScale', 1.20);

    assert(info.is_complete_in_display, ...
        'complete_molecule_ids returned a molecule not complete in display.');

    assert(info.n_fragments == 1, ...
        'Complete molecules should have one displayed fragment.');

    assert(abs(info.largest_fragment_fraction - 1.0) < 1e-12, ...
        'Complete molecules should have largest_fragment_fraction = 1.');
end

% There should also be at least one incomplete component at the displayed
% outer boundary for this wrapped template.
assert(any(~T.is_complete_in_display), ...
    'Boundary-crossing test should produce at least one incomplete boundary component.');

% Check that no incomplete molecule ID leaks into complete_molecule_ids.
incompleteIDs = T.molecule_id(~T.is_complete_in_display);

assert(isempty(intersect(completeIDs(:), incompleteIDs(:))), ...
    'complete_molecule_ids should not return incomplete components.');

% Provenance fields should still exist, but they are not canonical molecule
% IDs. A reconstructed molecule may span more than one cell_shift.
assert(numel(sys.base_mol_id) == sys.n_sites, ...
    'base_mol_id provenance should have one entry per site.');

assert(size(sys.cell_shift, 1) == sys.n_sites && size(sys.cell_shift, 2) == 3, ...
    'cell_shift provenance should be n_sites x 3.');

% At least one complete molecule should span more than one replicated cell
% shift along x. That demonstrates why [base_mol_id, cell_shift] alone is
% insufficient for molecule identity in this boundary-crossing case.
spansCellShift = false;

for k = 1:numel(completeIDs)
    idx = builder.site_indices_for_molecule(sys, completeIDs(k));
    shifts = unique(sys.cell_shift(idx, :), 'rows');

    if size(shifts, 1) > 1
        spansCellShift = true;
        break;
    end
end

assert(spansCellShift, ...
    ['At least one reconstructed complete molecule should span multiple ' ...
     'cell_shift values in the boundary-crossing case.']);

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

function [filename, tmpDir] = local_write_boundary_propene_poscar()
tmpDir = tempname;
mkdir(tmpDir);

filename = fullfile(tmpDir, 'POSCAR_boundary_propene');

% Cartesian Angstrom coordinates in a 20 A box.
%
% This is the same propene-like connectivity as the other builder tests, but
% shifted so that the C1 side of the molecule wraps across the x boundary:
%
%   C1 near x = 19.0
%   C2 near x = 0.45
%   C3 near x = 1.90
%
% Minimum-image C1-C2 distance is 1.45 A.
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

xyz = [
    19.00  5.00  5.00   % C1
     0.45  5.00  5.00   % C2
     1.90  5.00  5.00   % C3
    18.55  5.95  5.00   % H on C1
    18.55  4.05  5.00   % H on C1
     0.45  6.09  5.00   % H on C2
     2.35  5.95  5.00   % H on C3
     2.35  4.05  5.00   % H on C3
     1.90  5.00  6.09   % H on C3
];

fid = fopen(filename, 'w');
assert(fid > 0, 'Failed to open temporary POSCAR file.');

fprintf(fid, 'Boundary-crossing propene-like builder test\n');
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