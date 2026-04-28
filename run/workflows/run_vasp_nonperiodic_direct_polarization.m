%% run_vasp_nonperiodic_direct_polarization
% Real VASP nonperiodic direct-polarization workflow:
%
%   VASP/CONTCAR
%   -> crystal template
%   -> supercell system
%   -> centered same-stack pair
%   -> uniform +1/-1 molecular charges
%   -> nonperiodic external field
%   -> SCF problem
%   -> polarization operator factory
%   -> direct Thole SCF
%   -> active-space polarization energy
%   -> selected-pair + induced-dipole visualization
%
% This is intentionally the dense/direct finite-cluster baseline. No caches,
% no MEX, no iterative solver, no periodic operator.

clear; clc; close all;

fprintf('\n============================================================\n');
fprintf('Real VASP nonperiodic direct polarization workflow\n');
fprintf('============================================================\n');

%% User controls

cfg = struct();

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

cfg.field = struct();
cfg.field.mode = 'nonperiodic';
cfg.field.exclude_self = true;
cfg.field.use_thole_damping = true;

cfg.scf = struct();
cfg.scf.tol = 1e-10;
cfg.scf.maxIter = 200;
cfg.scf.mixing = 1.0;
cfg.scf.omega = 1.0;
cfg.scf.verbose = true;

cfg.operator = struct();
cfg.operator.mode = 'nonperiodic';
cfg.operator.backend = 'dense';
cfg.operator.use_thole = true;
cfg.operator.softening = 0.0;
cfg.operator.rcut = 15.0;   % Inf = full dense baseline; 15 bohr is much faster
cfg.operator.verbose = true;

cfg.plot = struct();
cfg.plot.showCOM = true;
cfg.plot.showLabels = true;
cfg.plot.showPairLine = true;
cfg.plot.drawBox = true;

cfg.plot.showDipoles = true;
cfg.plot.onlyPolarizableDipoles = true;
cfg.plot.maxArrows = 300;
cfg.plot.dipoleThreshold = 0.0;
cfg.plot.arrowScale = 120.0;   % visual scale only
cfg.plot.arrowLineWidth = 1.0;

HARTREE_TO_EV = 27.211386245988;

%% Check file

if ~isfile(cfg.filename)
    error('VASP file not found:\n  %s\nEdit cfg.filename at the top of this script.', cfg.filename);
end

fprintf('\nInput structure:\n  %s\n', cfg.filename);

%% 1. Import crystal template

fprintf('\n[1] Importing crystal template...\n');

crystal = io.import_contcar_as_crystal(cfg.filename, ...
    'BondScale', cfg.bondScale, ...
    'SortMolecules', false);

fprintf('  nSites     = %d\n', crystal.nSites);
fprintf('  nBaseMols  = %d\n', crystal.nBaseMols);
fprintf('  units      = %s\n', crystal.units.length);

%% 2. Model

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

%% 3. Build system

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

%% 4. Select centered neighbor pair

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

%% 5. Apply uniform charges

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

%% 6. Extract nonperiodic polarization system

fprintf('\n[6] Extracting nonperiodic polarization system...\n');

polsys = builder.extract_polarization_system(sys, struct('mode', 'nonperiodic'));
io.assert_atomic_units(polsys);

fprintf('  polsys.n_sites       = %d\n', polsys.n_sites);
fprintf('  polsys.is_periodic   = %d\n', polsys.is_periodic);
fprintf('  polarizable sites    = %d\n', nnz(polsys.site_is_polarizable));
fprintf('  charged sites        = %d\n', nnz(abs(polsys.site_charge) > 0));

%% 7. Compute nonperiodic external field

fprintf('\n[7] Computing nonperiodic external field...\n');

fieldParams = struct();
fieldParams.use_thole = cfg.field.use_thole_damping;
fieldParams.field = struct();
fieldParams.field.mode = cfg.field.mode;
fieldParams.field.exclude_self = cfg.field.exclude_self;
fieldParams.field.use_thole_damping = cfg.field.use_thole_damping;
fieldParams.field.target_mask = logical(polsys.site_is_polarizable(:));
fieldParams.field.source_mask = abs(polsys.site_charge(:)) > 0;

tField = tic;
Eext = calc.compute_external_field(polsys, fieldParams);
fieldTime = toc(tField);

fprintf('  Eext computed in %.6f s\n', fieldTime);
fprintf('  ||Eext||_F = %.12e\n', norm(Eext, 'fro'));
fprintf('  ||Eext polarizable||_F = %.12e\n', norm(Eext(polsys.site_is_polarizable, :), 'fro'));

%% 8. Prepare SCF problem

fprintf('\n[8] Preparing SCF problem...\n');

problem = thole.prepare_scf_problem(polsys, Eext, cfg.scf);

fprintf('  nPolSites = %d\n', problem.nPolSites);
fprintf('  active vector length = %d\n', numel(problem.Eext_pol_vec));

%% 9. Build polarization operator

fprintf('\n[9] Building polarization operator...\n');

tOp = tic;
op = thole.make_polarization_operator(polsys, problem, ...
    'Mode', cfg.operator.mode, ...
    'Backend', cfg.operator.backend, ...
    'UseThole', cfg.operator.use_thole, ...
    'Softening', cfg.operator.softening, ...
    'Rcut', cfg.operator.rcut, ...
    'Verbose', cfg.operator.verbose);
opTime = toc(tOp);

fprintf('  operator build wrapper time = %.6f s\n', opTime);
fprintf('  operator mode               = %s\n', op.mode);
fprintf('  operator backend            = %s\n', op.backend);
fprintf('  operator kind               = %s\n', op.kind);
fprintf('  operator size               = %d x %d\n', op.size(1), op.size(2));

if isfield(op, 'info')
    fprintf('  opinfo assembly time        = %.6f s\n', op.info.assembly_time);
    fprintf('  pair blocks total           = %d\n', op.info.nPairBlocks);
    fprintf('  pair blocks kept            = %d\n', op.info.nPairBlocksKept);
    fprintf('  skipped by cutoff           = %d\n', op.info.nPairBlocksSkippedCutoff);
    fprintf('  rcut                        = %.6g bohr\n', op.info.rcut);
end

%% 10. Direct SCF solve

fprintf('\n[10] Solving direct SCF...\n');

solveOpts = struct();
solveOpts.compute_residual = true;

tSolve = tic;
[mu, solveInfo] = thole.solve_scf_direct(problem, op, solveOpts);
solveTime = toc(tSolve);

fprintf('  solve call time   = %.6f s\n', solveTime);
fprintf('  solver total time = %.6f s\n', solveInfo.total_time);
fprintf('  setup time        = %.6f s\n', solveInfo.setup_time);
fprintf('  linear solve time = %.6f s\n', solveInfo.solve_time);
fprintf('  residual enabled  = %d\n', solveInfo.compute_residual);

if solveInfo.compute_residual
    fprintf('  residual time     = %.6f s\n', solveInfo.residual_time);
    fprintf('  info.relres       = %.12e\n', solveInfo.relres);
else
    fprintf('  info.relres       = skipped\n');
end

fprintf('  ||mu||_F = %.12e\n', norm(mu, 'fro'));
fprintf('  max |mu_i| = %.12e\n', max(vecnorm(mu, 2, 2)));

%% 11. Energy

fprintf('\n[11] Computing active-space polarization energy...\n');

energy = calc.compute_total_energy_active_space(polsys, problem, mu, Eext, op);

fprintf('\nEnergy breakdown:\n');
fprintf('  polarization_self       = %+ .12e Ha  (%+ .8f eV)\n', ...
    energy.polarization_self, energy.polarization_self * HARTREE_TO_EV);
fprintf('  external_charge_dipole  = %+ .12e Ha  (%+ .8f eV)\n', ...
    energy.external_charge_dipole, energy.external_charge_dipole * HARTREE_TO_EV);
fprintf('  dipole_dipole           = %+ .12e Ha  (%+ .8f eV)\n', ...
    energy.dipole_dipole, energy.dipole_dipole * HARTREE_TO_EV);
fprintf('  total                   = %+ .12e Ha  (%+ .8f eV)\n', ...
    energy.total, energy.total * HARTREE_TO_EV);
fprintf('  total_stationary        = %+ .12e Ha  (%+ .8f eV)\n', ...
    energy.total_stationary, energy.total_stationary * HARTREE_TO_EV);
fprintf('  stationary_consistency  = %+ .12e Ha\n', energy.stationary_consistency);

if abs(energy.stationary_consistency) > 1e-8
    warning('Stationary energy consistency is larger than expected.');
end

%% 12. Visualize selected pair and induced dipoles

fprintf('\n[12] Plotting selected pair and induced dipoles...\n');

fig = figure('Name', 'Nonperiodic direct polarization');
ax = axes(fig);

plotTitle = sprintf('Nonperiodic direct: E_{pol} = %.4f eV', energy.total * HARTREE_TO_EV);

plotOut = viz.plot_supercell_selection(sys, selection, ...
    'Axes', ax, ...
    'Title', plotTitle, ...
    'DrawBox', cfg.plot.drawBox, ...
    'ShowCOM', cfg.plot.showCOM, ...
    'ShowLabels', cfg.plot.showLabels, ...
    'ShowPairLine', cfg.plot.showPairLine);

if cfg.plot.showDipoles
    arrowOut = local_plot_induced_dipoles(ax, polsys, mu, cfg.plot);

    fprintf('\nDipole visualization:\n');
    fprintf('  arrows plotted = %d\n', arrowOut.n_arrows);
    fprintf('  max |mu|       = %.12e\n', arrowOut.max_mu_norm);
end

fprintf('\nVisualization output:\n');
fprintf('  ref molecule ID      = %d\n', plotOut.ref_mol_id);
fprintf('  neighbor molecule ID = %d\n', plotOut.neighbor_mol_id);
fprintf('  pair distance / bohr = %.6f\n', plotOut.pair_distance);
fprintf('  pair midpoint / bohr:\n');
disp(plotOut.pair_midpoint);

fprintf('\nWorkflow completed successfully.\n');

%% Local helpers

function out = local_plot_induced_dipoles(ax, polsys, mu, plotCfg)

muNorm = vecnorm(mu, 2, 2);

if plotCfg.onlyPolarizableDipoles
    mask = logical(polsys.site_is_polarizable(:));
else
    mask = true(polsys.n_sites, 1);
end

mask = mask & (muNorm(:) > plotCfg.dipoleThreshold);

idx = find(mask);

if isempty(idx)
    out = struct();
    out.handle = gobjects(0);
    out.indices = [];
    out.n_arrows = 0;
    out.max_mu_norm = 0;
    return;
end

[~, order] = sort(muNorm(idx), 'descend');
idx = idx(order);

if numel(idx) > plotCfg.maxArrows
    idx = idx(1:plotCfg.maxArrows);
end

X = polsys.site_pos(idx, :);
U = mu(idx, :) * plotCfg.arrowScale;

hold(ax, 'on');

h = quiver3(ax, ...
    X(:,1), X(:,2), X(:,3), ...
    U(:,1), U(:,2), U(:,3), ...
    0, ...
    'k', ...
    'LineWidth', plotCfg.arrowLineWidth, ...
    'MaxHeadSize', 0.6);

out = struct();
out.handle = h;
out.indices = idx(:);
out.n_arrows = numel(idx);
out.max_mu_norm = max(muNorm(idx));

end