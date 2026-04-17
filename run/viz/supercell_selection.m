%% run_test_plot_supercell_selection
clear; clc;
close all;

rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
filename   = fullfile(rootFolder, 'a_0.0_CONTCAR.vasp');

%% ------------------------------------------------------------------------
% Import crystal template
crystal = io.import_contcar_as_crystal(filename, ...
    'BondScale', 1.20, ...
    'WrapFractional', true, ...
    'SortMolecules', false);

fprintf('Imported crystal template:\n');
fprintf('  nSites     = %d\n', size(crystal.cart_coords, 1));
fprintf('  nBaseMols  = %d\n', numel(unique(crystal.base_mol_id)));

%% ------------------------------------------------------------------------
% Minimal model
model = struct();
model.polarizable_types = {'H','C','N','O'};
model.alpha_by_type = struct( ...
    'H', 0.696, ...
    'C', 1.750, ...
    'N', 1.073, ...
    'O', 0.837);
model.thole_a = 0.39;

%% ------------------------------------------------------------------------
% Working options
opts = struct();
opts.supercell_size = [2 5 1];
opts.bondScale = 1.20;
opts.verbose = true;
opts.removeMolIDs = [];
opts.activeMolIDs = [];

%% ------------------------------------------------------------------------
% Build working system
sys = builder.make_crystal_system(crystal, model, opts);

fprintf('\nBuilt working system:\n');
fprintf('  n_unit_sites = %d\n', sys.n_unit_sites);
fprintf('  n_cells      = %d\n', sys.n_cells);
fprintf('  n_sites      = %d\n', sys.n_sites);

ccompleteIDs = builder.complete_molecule_ids(sys, 'Verbose', true);

[refMolID, refSummary] = builder.choose_center_reference_molecule(sys, ...
    'RequireComplete', true, ...
    'Verbose', true);

desc = builder.complete_molecule_descriptors_relative_to_reference(sys, refMolID, ...
    'StackAxis', 'b', ...
    'Verbose', true);

disp(desc.table(1:min(12,height(desc.table)), :))

%% ------------------------------------------------------------------------
% Plot
fig = figure(1);
clf(fig);
fig.Color = 'w';
ax = axes('Parent', fig);

out = viz.plot_supercell_selection(sys, ...
    'ReferenceMolID', refMolID, ...
    'RequireCompleteRef', true, ...
    'IncompleteAction', 'error', ...
    'EnvMarkerSize', 52, ...
    'RefMarkerSize', 78, ...
    'ShowCOM', true, ...
    'ShowLabels', true, ...
    'DrawBox', true, ...
    'Axes', ax, ...
    'Title', 'Supercell with centered complete reference molecule');

if out.reference_unique_mol_id ~= refMolID
    error('Plotter reference ID mismatch.');
end

fprintf('\nPlot complete.\n');
fprintf('  plotted reference ID = %d\n', out.reference_unique_mol_id);

fprintf('\nCheckpoint passed.\n');