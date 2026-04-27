%% tutorial_04_crystal_builder_neighbor_selection_and_viz
% Demonstrates:
%   1. Build a toy VASP crystal
%   2. Import as a crystal template
%   3. Build a molecule-aware supercell
%   4. Choose centered neighbor pairs by relation
%   5. Plot selected pairs
%
% Toy packing:
%   - lower slab: 4 propene molecules in a 2x2 motif
%   - upper slab: same motif shifted upward in z
%
% This gives:
%   - same-stack neighbors (along +b / y)
%   - side-stack family 1
%   - side-stack family 2
%
% Assumes src/ is already on the MATLAB path.

clear; clc;

fprintf('\n============================================================\n');
fprintf('Tutorial 04: crystal builder, neighbor selection, and viz\n');
fprintf('============================================================\n');

%% 1. Write a temporary POSCAR-style file

tmpDir = tempname;
mkdir(tmpDir);
cleanup = onCleanup(@() local_cleanup(tmpDir));

filename = fullfile(tmpDir, 'POSCAR_propene_builder_viz');

% Approximate propene geometry in Angstrom.
%
% Connectivity:
%   C1H2 = C2H - C3H3
%
% Approximate geometry:
%   C1=C2       ~1.34 A
%   C2-C3       ~1.50 A
%   C1-C2-C3    ~120 deg
%   C-H         ~1.09 A
%
% Methyl hydrogens are tetrahedral around C3 and oriented away from C2.
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
xyz0 = [
    5.0000  5.0000  5.0000   % C1
    6.3400  5.0000  5.0000   % C2
    7.0900  6.2990  5.0000   % C3

    4.4550  5.9439  5.0000   % H on C1
    4.4550  4.0561  5.0000   % H on C1
    6.8850  4.0561  5.0000   % H on C2

    6.3817  7.1275  5.0000   % H on C3
    7.7166  6.3568  5.8900   % H on C3
    7.7166  6.3568  4.1100   % H on C3
];

% Packing shifts in Angstrom.
shiftSameStack = [0.0  6.0  0.0];   % along b / y
shiftSideStack = [7.0  0.0  0.0];   % along a / x
shiftUpperSlab = [0.0  0.0  4.0];   % upper slab in z

% Lower slab
xyzA = xyz0;
xyzB = xyz0 + shiftSameStack;
xyzC = xyz0 + shiftSideStack;
xyzD = xyz0 + shiftSideStack + shiftSameStack;

% Upper slab
xyzE = xyzA + shiftUpperSlab;
xyzF = xyzB + shiftUpperSlab;
xyzG = xyzC + shiftUpperSlab;
xyzH = xyzD + shiftUpperSlab;

allMols = {xyzA, xyzB, xyzC, xyzD, xyzE, xyzF, xyzG, xyzH};

% VASP requires species-grouped coordinates.
xyzC_all = [];
xyzH_all = [];
for m = 1:numel(allMols)
    xyzM = allMols{m};
    xyzC_all = [xyzC_all; xyzM(1:3, :)]; %#ok<AGROW>
    xyzH_all = [xyzH_all; xyzM(4:9, :)]; %#ok<AGROW>
end
xyz = [xyzC_all; xyzH_all];

fid = fopen(filename, 'w');
if fid < 0
    error('Failed to open temporary POSCAR file.');
end

fprintf(fid, 'Eight propene-like molecules for builder/viz tutorial\n');
fprintf(fid, '1.0\n');

% A somewhat compact, more cube-like unit cell.
fprintf(fid, '15.0 0.0 0.0\n');
fprintf(fid, '0.0 12.0 0.0\n');
fprintf(fid, '0.0 0.0 10.0\n');

fprintf(fid, 'C H\n');
fprintf(fid, '24 48\n');
fprintf(fid, 'Cartesian\n');

for i = 1:size(xyz, 1)
    fprintf(fid, '%.10f %.10f %.10f\n', xyz(i,1), xyz(i,2), xyz(i,3));
end
fclose(fid);

fprintf('\nWrote temporary POSCAR:\n  %s\n', filename);

%% 2. Import as a crystal template

crystal = io.import_contcar_as_crystal(filename, ...
    'BondScale', 1.20, ...
    'SortMolecules', false);

fprintf('\nCrystal template:\n');
fprintf('  nSites       = %d\n', crystal.nSites);
fprintf('  nBaseMols    = %d\n', crystal.nBaseMols);
fprintf('  units.length = %s\n', crystal.units.length);

fprintf('\nBase molecule IDs in unit cell:\n');
disp(unique(crystal.base_mol_id(:)).');

%% 3. Define a simple polarizability model

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
model.alpha_by_class.C_deg3 = 1.334;
model.alpha_by_class.C_deg4 = 1.750;
model.alpha_by_class.H_on_C_deg3 = 0.496;
model.alpha_by_class.H_on_C_deg4 = 0.696;

fprintf('\nPolarizability model:\n');
fprintf('  alpha input units = %s\n', model.alpha_units);
fprintf('  C alpha input     = %.3f Angstrom^3\n', model.alpha_by_class.C_deg3);
fprintf('  H alpha input     = %.3f Angstrom^3\n', model.alpha_by_class.H_on_C_deg3);

%% 4. Build a supercell

reps = [2 2 2];

sys = builder.make_crystal_system(crystal, model, struct( ...
    'supercell_size', reps, ...
    'bondScale', 1.20, ...
    'verbose', false));

fprintf('\nBuilt crystal system:\n');
fprintf('  supercell_size = [%d %d %d]\n', sys.supercell_size);
fprintf('  n_sites        = %d\n', sys.n_sites);
fprintf('  n molecules    = %d\n', numel(sys.molecule_table.molecule_id));
fprintf('  n complete     = %d\n', nnz(sys.molecule_table.is_complete_in_display));
fprintf('  length unit    = %s\n', sys.units.length);
fprintf('  alpha unit     = %s\n', sys.units.alpha);

completeIDs = builder.complete_molecule_ids(sys);

fprintf('\nComplete molecule IDs:\n');
disp(completeIDs(:).');

%% 5. Choose a centered reference molecule

[refMolID, refSummary] = builder.choose_center_reference_molecule(sys, ...
    'RequireComplete', true, ...
    'Verbose', false);

fprintf('\nCentered reference molecule:\n');
fprintf('  refMolID = %d\n', refMolID);
fprintf('  supercell center / bohr:\n');
disp(refSummary.supercell_center);
fprintf('  reference COM / bohr:\n');
disp(refSummary.chosen_com);

%% 6. Build descriptors relative to the reference

desc = builder.complete_molecule_descriptors_relative_to_reference(sys, refMolID, ...
    'StackAxis', 'b', ...
    'IncludeReference', false, ...
    'IncludeNormals', true);

fprintf('\nDescriptor table relative to reference %d:\n', refMolID);
disp(desc.table(:, {'molecule_id', 'd_par', 'd_perp', 'distance', 'normal_angle_deg'}));

%% 7. Select a same-stack pair

sameStackSelection = builder.select_centered_neighbor_pair(sys, ...
    'Relation', 'same_stack', ...
    'Shell', 1, ...
    'StackAxis', 'b', ...
    'Direction', 'either', ...
    'PerpTol', 1e-6, ...
    'ShellTol', 1e-6, ...
    'Verbose', false);

fprintf('\nSelected same-stack pair:\n');
fprintf('  relation       = %s\n', sameStackSelection.relation);
fprintf('  reference ID   = %d\n', sameStackSelection.reference_mol_id);
fprintf('  neighbor ID    = %d\n', sameStackSelection.neighbor_mol_id);
fprintf('  midpoint distance from supercell center / bohr = %.6f\n', sameStackSelection.midpoint_distance);
fprintf('  pair midpoint / bohr:\n');
disp(sameStackSelection.pair_midpoint);

%% 8. Select two side-stack families

sideStackSelection1 = builder.select_centered_neighbor_pair(sys, ...
    'Relation', 'side_stack', ...
    'Shell', 1, ...
    'Member', 1, ...
    'StackAxis', 'b', ...
    'Direction', 'either', ...
    'SameStackPerpTol', 1e-6, ...
    'ShellTol', 1e-6, ...
    'VectorGroupTol', 1e-6, ...
    'Verbose', false);

fprintf('\nSelected side-stack family 1:\n');
fprintf('  relation       = %s\n', sideStackSelection1.relation);
fprintf('  reference ID   = %d\n', sideStackSelection1.reference_mol_id);
fprintf('  neighbor ID    = %d\n', sideStackSelection1.neighbor_mol_id);
fprintf('  midpoint distance from supercell center / bohr = %.6f\n', sideStackSelection1.midpoint_distance);
fprintf('  pair midpoint / bohr:\n');
disp(sideStackSelection1.pair_midpoint);

sideStackSelection2 = builder.select_centered_neighbor_pair(sys, ...
    'Relation', 'side_stack', ...
    'Shell', 1, ...
    'Member', 2, ...
    'StackAxis', 'b', ...
    'Direction', 'either', ...
    'SameStackPerpTol', 1e-6, ...
    'ShellTol', 1e-6, ...
    'VectorGroupTol', 1e-6, ...
    'Verbose', false);

fprintf('\nSelected side-stack family 2:\n');
fprintf('  relation       = %s\n', sideStackSelection2.relation);
fprintf('  reference ID   = %d\n', sideStackSelection2.reference_mol_id);
fprintf('  neighbor ID    = %d\n', sideStackSelection2.neighbor_mol_id);
fprintf('  midpoint distance from supercell center / bohr = %.6f\n', sideStackSelection2.midpoint_distance);
fprintf('  pair midpoint / bohr:\n');
disp(sideStackSelection2.pair_midpoint);

%% 9. Visualize the three selected pairs

fprintf('\nPlotting selected pairs...\n');

fig1 = figure('Name', 'Tutorial 04: same-stack pair');
ax1 = axes(fig1);

samePlot = viz.plot_supercell_selection(sys, sameStackSelection, ...
    'Axes', ax1, ...
    'Title', 'Selected same-stack pair', ...
    'DrawBox', true, ...
    'ShowCOM', true, ...
    'ShowLabels', true, ...
    'ShowPairLine', true);

fprintf('\nSame-stack visualization output:\n');
fprintf('  ref molecule ID      = %d\n', samePlot.ref_mol_id);
fprintf('  neighbor molecule ID = %d\n', samePlot.neighbor_mol_id);
fprintf('  pair distance / bohr = %.6f\n', samePlot.pair_distance);
fprintf('  pair midpoint / bohr:\n');
disp(samePlot.pair_midpoint);

fig2 = figure('Name', 'Tutorial 04: side-stack family 1');
ax2 = axes(fig2);

sidePlot1 = viz.plot_supercell_selection(sys, sideStackSelection1, ...
    'Axes', ax2, ...
    'Title', 'Selected side-stack family 1', ...
    'DrawBox', true, ...
    'ShowCOM', true, ...
    'ShowLabels', true, ...
    'ShowPairLine', true);

fprintf('\nSide-stack family 1 visualization output:\n');
fprintf('  ref molecule ID      = %d\n', sidePlot1.ref_mol_id);
fprintf('  neighbor molecule ID = %d\n', sidePlot1.neighbor_mol_id);
fprintf('  pair distance / bohr = %.6f\n', sidePlot1.pair_distance);
fprintf('  pair midpoint / bohr:\n');
disp(sidePlot1.pair_midpoint);

fig3 = figure('Name', 'Tutorial 04: side-stack family 2');
ax3 = axes(fig3);

sidePlot2 = viz.plot_supercell_selection(sys, sideStackSelection2, ...
    'Axes', ax3, ...
    'Title', 'Selected side-stack family 2', ...
    'DrawBox', true, ...
    'ShowCOM', true, ...
    'ShowLabels', true, ...
    'ShowPairLine', true);

fprintf('\nSide-stack family 2 visualization output:\n');
fprintf('  ref molecule ID      = %d\n', sidePlot2.ref_mol_id);
fprintf('  neighbor molecule ID = %d\n', sidePlot2.neighbor_mol_id);
fprintf('  pair distance / bohr = %.6f\n', sidePlot2.pair_distance);
fprintf('  pair midpoint / bohr:\n');
disp(sidePlot2.pair_midpoint);

fprintf('\nTutorial 04 completed successfully.\n');

function local_cleanup(tmpDir)
if exist(tmpDir, 'dir')
    rmdir(tmpDir, 's');
end
end