%% run_real_vasp_uniform_charged_pair_viz
% Real-structure smoke run:
%
%   VASP/CONTCAR
%   -> crystal template
%   -> supercell system
%   -> centered same-stack pair
%   -> uniform +1/-1 molecular charges
%   -> charged pair removed from polarizable mask
%   -> selected charged pair visualized
%
% This is meant to validate that the cleaned builder/viz stack reproduces
% the visual style of the previous induced-dipole comparison script:
%   grey environment, bold selected pair, COMs, labels, pair line, box.

clear; clc; close all;

fprintf('\n============================================================\n');
fprintf('Real VASP uniform-charged pair visualization\n');
fprintf('============================================================\n');

%% User controls

cfg = struct();

% Same default location used in the earlier induced-dipole visualization run.
cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

cfg.supercellSize = [2 5 1];
cfg.bondScale = 1.20;

cfg.relation = 'same_stack';
cfg.shell = 1;
cfg.stackAxis = 'b';
cfg.direction = 'either';

cfg.pairCharges = [+1 -1];

cfg.verbose = true;

% Plot style: intentionally close to the older induced-dipole visual style.
cfg.plot.showCOM = true;
cfg.plot.showLabels = true;
cfg.plot.showPairLine = true;
cfg.plot.drawBox = true;

%% Check file

if ~isfile(cfg.filename)
    error('VASP file not found:\n  %s\nEdit cfg.filename at the top of this script.', cfg.filename);
end

fprintf('\nInput structure:\n  %s\n', cfg.filename);

%% Import crystal template

fprintf('\n[1] Importing crystal template...\n');

crystal = io.import_contcar_as_crystal(cfg.filename, ...
    'BondScale', cfg.bondScale, ...
    'SortMolecules', false);

fprintf('  nSites     = %d\n', crystal.nSites);
fprintf('  nBaseMols  = %d\n', crystal.nBaseMols);
fprintf('  units      = %s\n', crystal.units.length);

%% Model

fprintf('\n[2] Defining polarizability model...\n');

model = struct();

model.thole_a = 0.39;
model.alpha_units = 'angstrom^3';

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

fprintf('  alpha_units = %s\n', model.alpha_units);
fprintf('  thole_a     = %.3f\n', model.thole_a);

%% Build system

fprintf('\n[3] Building crystal system...\n');

buildOpts = struct();
buildOpts.supercell_size = cfg.supercellSize;
buildOpts.bondScale = cfg.bondScale;
buildOpts.verbose = cfg.verbose;

sys0 = builder.make_crystal_system(crystal, model, buildOpts);

fprintf('\nBuilt system:\n');
fprintf('  supercell_size = [%d %d %d]\n', sys0.supercell_size);
fprintf('  n_sites        = %d\n', sys0.n_sites);
fprintf('  n_molecules    = %d\n', numel(sys0.molecule_table.molecule_id));
fprintf('  n_complete     = %d\n', nnz(sys0.molecule_table.is_complete_in_display));
fprintf('  length unit    = %s\n', sys0.units.length);
fprintf('  alpha unit     = %s\n', sys0.units.alpha);

completeIDs = builder.complete_molecule_ids(sys0);

if isempty(completeIDs)
    error('No complete molecules found in the displayed supercell.');
end

fprintf('\nComplete molecule IDs:\n');
disp(completeIDs(:).');

%% Select centered same-stack pair

fprintf('\n[4] Selecting centered %s shell %d pair...\n', cfg.relation, cfg.shell);

selection = builder.select_centered_neighbor_pair(sys0, ...
    'Relation', cfg.relation, ...
    'Shell', cfg.shell, ...
    'StackAxis', cfg.stackAxis, ...
    'Direction', cfg.direction, ...
    'Verbose', cfg.verbose);

pairVector0 = selection.neighbor_com - selection.reference_com;
pairDistance0 = norm(pairVector0);

fprintf('\nSelected pair before charging:\n');
fprintf('  relation       = %s\n', selection.relation);
fprintf('  reference ID   = %d\n', selection.reference_mol_id);
fprintf('  neighbor ID    = %d\n', selection.neighbor_mol_id);
fprintf('  pair distance  = %.6f bohr\n', pairDistance0);
fprintf('  midpoint dist  = %.6f bohr\n', selection.midpoint_distance);
fprintf('  pair midpoint / bohr:\n');
disp(selection.pair_midpoint);

%% Apply uniform charges

fprintf('\n[5] Applying uniform charges...\n');

molIDs = [selection.reference_mol_id selection.neighbor_mol_id];

sys = builder.apply_molecule_charges(sys0, molIDs, ...
    'Mode', 'uniform', ...
    'TotalCharges', cfg.pairCharges, ...
    'SetActive', true, ...
    'DisablePolarizabilityOnCharged', true, ...
    'ZeroExistingCharges', true, ...
    'RequireComplete', true, ...
    'Verbose', cfg.verbose);

idxRef = builder.site_indices_for_molecule(sys, selection.reference_mol_id);
idxNbr = builder.site_indices_for_molecule(sys, selection.neighbor_mol_id);

fprintf('\nCharge summary:\n');
fprintf('  ref molecule %d total charge = %+ .8f e\n', ...
    selection.reference_mol_id, sum(sys.site_charge(idxRef)));
fprintf('  nbr molecule %d total charge = %+ .8f e\n', ...
    selection.neighbor_mol_id, sum(sys.site_charge(idxNbr)));
fprintf('  total system charge          = %+ .8f e\n', sum(sys.site_charge));

chargedMask = ismember(sys.site_mol_id, molIDs);

fprintf('\nActive/polarizable summary:\n');
fprintf('  active sites                  = %d\n', nnz(sys.site_is_active));
fprintf('  charged sites                 = %d\n', nnz(chargedMask));
fprintf('  charged sites polarizable?    = %d\n', any(sys.site_is_polarizable(chargedMask)));
fprintf('  remaining polarizable sites   = %d\n', nnz(sys.site_is_polarizable));

if any(sys.site_is_polarizable(chargedMask))
    error('Charged molecule sites should not remain polarizable.');
end

io.assert_atomic_units(sys);

%% Extract polsys as a quick solver-facing check

fprintf('\n[6] Extracting periodic polarization system...\n');

polsys = builder.extract_polarization_system(sys, struct('mode', 'periodic'));
io.assert_atomic_units(polsys);

fprintf('  polsys.n_sites       = %d\n', polsys.n_sites);
fprintf('  polsys.is_periodic   = %d\n', polsys.is_periodic);
fprintf('  polarizable sites    = %d\n', nnz(polsys.site_is_polarizable));
fprintf('  charged sites        = %d\n', nnz(abs(polsys.site_charge) > 0));

%% Visualize charged selected pair

fprintf('\n[7] Plotting charged selected pair...\n');

fig = figure('Name', 'Real VASP uniform-charged pair');
ax = axes(fig);

plotOut = viz.plot_supercell_selection(sys, selection, ...
    'Axes', ax, ...
    'Title', 'Uniform charged same-stack pair', ...
    'DrawBox', cfg.plot.drawBox, ...
    'ShowCOM', cfg.plot.showCOM, ...
    'ShowLabels', cfg.plot.showLabels, ...
    'ShowPairLine', cfg.plot.showPairLine);

fprintf('\nVisualization output:\n');
fprintf('  ref molecule ID      = %d\n', plotOut.ref_mol_id);
fprintf('  neighbor molecule ID = %d\n', plotOut.neighbor_mol_id);
fprintf('  pair distance / bohr = %.6f\n', plotOut.pair_distance);
fprintf('  pair midpoint / bohr:\n');
disp(plotOut.pair_midpoint);

fprintf('\nRun completed successfully.\n');