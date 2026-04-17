%% run_test_build_uniform_charged_pair_system
clear; clc; close all;

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
% Model with simple graph-based polarizability classes
model = struct();

model.polarizable_classes = { ...
    'H_on_C_deg3', ...
    'H_on_C_deg4', ...
    'C_deg3', ...
    'C_deg4', ...
    'N', ...
    'O'};

model.alpha_by_class = struct( ...
    'H_on_C_deg3', 0.496, ...
    'H_on_C_deg4', 0.696, ...
    'C_deg3',      1.334, ...
    'C_deg4',      1.750, ...
    'N',           1.073, ...
    'O',           0.837);

model.thole_a = 0.39;

%% ------------------------------------------------------------------------
% Build working supercell
opts = struct();
opts.supercell_size = [2 5 1];
opts.bondScale = 1.20;
opts.verbose = true;
opts.removeMolIDs = [];
opts.activeMolIDs = [];

sys = builder.make_crystal_system(crystal, model, opts);

fprintf('\nBuilt working system:\n');
fprintf('  n_unit_sites = %d\n', sys.n_unit_sites);
fprintf('  n_cells      = %d\n', sys.n_cells);
fprintf('  n_sites      = %d\n', sys.n_sites);

%% ------------------------------------------------------------------------
% Example pair you already identified for this supercell:
% first-shell same-stack neighbor
refID = 32;
nbrID = 38;

fprintf('\nChosen pair:\n');
fprintf('  ref ID = %d\n', refID);
fprintf('  nbr ID = %d\n', nbrID);

%% ------------------------------------------------------------------------
% Apply uniform cation/anion charges and disable polarizability on those sites
sys = builder.apply_molecule_charges(sys, [refID nbrID], ...
    'Mode', 'uniform', ...
    'TotalCharges', [+1 -1], ...
    'SetActive', true, ...
    'DisablePolarizabilityOnCharged', true, ...
    'ZeroExistingCharges', true, ...
    'Verbose', true);

%% ------------------------------------------------------------------------
% Diagnostics
refIdx = builder.site_indices_for_molecule(sys, refID);
nbrIdx = builder.site_indices_for_molecule(sys, nbrID);

fprintf('\nPost-assignment diagnostics:\n');
fprintf('  ref site count                  = %d\n', numel(refIdx));
fprintf('  nbr site count                  = %d\n', numel(nbrIdx));
fprintf('  total system charge             = %+0.10f\n', sum(sys.site_charge));
fprintf('  ref total charge                = %+0.10f\n', sum(sys.site_charge(refIdx)));
fprintf('  nbr total charge                = %+0.10f\n', sum(sys.site_charge(nbrIdx)));
fprintf('  ref polarizable sites remaining = %d\n', nnz(sys.site_is_polarizable(refIdx)));
fprintf('  nbr polarizable sites remaining = %d\n', nnz(sys.site_is_polarizable(nbrIdx)));
fprintf('  active sites total              = %d\n', nnz(sys.site_is_active));

if abs(sum(sys.site_charge(refIdx)) - 1.0) > 1e-12
    error('Reference molecule charge assignment failed.');
end

if abs(sum(sys.site_charge(nbrIdx)) + 1.0) > 1e-12
    error('Neighbor molecule charge assignment failed.');
end

if any(sys.site_is_polarizable(refIdx))
    error('Reference molecule still has polarizable sites enabled.');
end

if any(sys.site_is_polarizable(nbrIdx))
    error('Neighbor molecule still has polarizable sites enabled.');
end

if ~all(sys.site_is_active(refIdx)) || ~all(sys.site_is_active(nbrIdx))
    error('Active mask was not set correctly on charged molecules.');
end

%% ------------------------------------------------------------------------
% Optional plot
fig = figure(1);
clf(fig);
fig.Color = 'w';
ax = axes('Parent', fig);

viz.plot_supercell_selection(sys, ...
    'ReferenceMolID', refID, ...
    'NeighborMolID', nbrID, ...
    'RequireCompleteRef', true, ...
    'RequireCompleteNbr', true, ...
    'IncompleteAction', 'error', ...
    'EnvMarkerSize', 52, ...
    'RefMarkerSize', 78, ...
    'NbrMarkerSize', 78, ...
    'ShowCOM', true, ...
    'ShowLabels', true, ...
    'DrawBox', true, ...
    'Axes', ax, ...
    'Title', sprintf('Uniform charged pair: ref %d (+1), nbr %d (-1)', refID, nbrID));

%% ------------------------------------------------------------------------
% Hand off to calculator
polsys = sys;

fprintf('\nCheckpoint passed.\n');
fprintf('polsys is ready to hand to the polarization calculator.\n');