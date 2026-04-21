%% compare_induced_dipoles_nonperiodic_vs_periodic
% Visualize induced dipoles for the same periodic polarization operator
% under:
%   (1) nonperiodic external field
%   (2) periodic external field
%
% Assumes your current Polarize refactor package layout.

clear; clc; close all;

%% ------------------------------------------------------------------------
% User controls

cfg = struct();

cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename   = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

cfg.supercellSize = [4 5 2];
cfg.bondScale     = 1.20;
cfg.pairCharges   = [+1 -1];

cfg.use_thole = true;
cfg.verbose   = true;
cfg.softening = 0.0;

cfg.periodic = struct();
cfg.periodic.alpha    = 0.30;
cfg.periodic.rcut     = 12.0;
cfg.periodic.kcut     = 3.5;
cfg.periodic.boundary = 'tinfoil';

cfg.sor = struct();
cfg.sor.tol           = 1e-6;
cfg.sor.maxIter       = 200;
cfg.sor.omega         = 0.45;
cfg.sor.printEvery    = 10;
cfg.sor.residualEvery = 10;
cfg.sor.stopMetric    = 'max_dmu';

cfg.kspace_mode = 'auto';   % use 'auto' if you want
cfg.k_block_size = 2048;
cfg.kspace_memory_limit_gb = 8;

% plotting controls
cfg.plot.onlyPolarizable = true;
cfg.plot.arrowScale = 100.0;      % tune visually
cfg.plot.maxArrows = 300;        % plot only largest dipoles
cfg.plot.showAllAtoms = true;
cfg.plot.atomMarkerSize = 8;
cfg.plot.pairMarkerSize = 40;
cfg.plot.dipoleThreshold = 0.0;  % can raise if cluttered

%% ------------------------------------------------------------------------
% Import crystal template

crystal = io.import_contcar_as_crystal(cfg.filename, ...
    'BondScale', cfg.bondScale, ...
    'WrapFractional', true, ...
    'SortMolecules', false);

%% ------------------------------------------------------------------------
% Model

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
% Build working system

buildOpts = struct();
buildOpts.supercell_size = cfg.supercellSize;
buildOpts.bondScale      = cfg.bondScale;
buildOpts.verbose        = cfg.verbose;
buildOpts.removeMolIDs   = [];
buildOpts.activeMolIDs   = [];

sys0 = builder.make_crystal_system(crystal, model, buildOpts);

%% ------------------------------------------------------------------------
% Choose charged pair automatically

[refID, ~] = builder.choose_center_reference_molecule(sys0, ...
    'RequireComplete', true, ...
    'Verbose', cfg.verbose);

desc = builder.complete_molecule_descriptors_relative_to_reference(sys0, refID, ...
    'Verbose', cfg.verbose);

sameStack = builder.select_same_stack_neighbor(desc, 1, ...
    'Direction', 'either', ...
    'Verbose', cfg.verbose);

if isempty(sameStack.match_table) || height(sameStack.match_table) < 1
    error('Automatic same-stack shell-1 neighbor selection failed.');
end

nbrID = sameStack.match_table.molecule_id(1);

fprintf('\nChosen pair:\n');
fprintf('  ref ID = %d\n', refID);
fprintf('  nbr ID = %d\n', nbrID);

%% ------------------------------------------------------------------------
% Apply charges and extract polarization system

sys = builder.apply_molecule_charges(sys0, [refID nbrID], ...
    'Mode', 'uniform', ...
    'TotalCharges', cfg.pairCharges, ...
    'SetActive', true, ...
    'DisablePolarizabilityOnCharged', true, ...
    'ZeroExistingCharges', true, ...
    'Verbose', cfg.verbose);

params_extract = struct();
params_extract.ewald = struct();
params_extract.ewald.mode = 'periodic';

polsys0 = builder.extract_polarization_system(sys, params_extract);
io.assert_atomic_units(polsys0);

H0 = local_get_direct_lattice(polsys0);
polsys = polsys0;
polsys.super_lattice = H0;

center0 = 0.5 * sum(H0, 2);
center_per = 0.5 * sum(polsys.super_lattice, 2);
shift = (center_per - center0).';
polsys.site_pos = polsys0.site_pos + shift;

ewaldParams = struct();
ewaldParams.alpha    = cfg.periodic.alpha;
ewaldParams.rcut     = cfg.periodic.rcut;
ewaldParams.kcut     = cfg.periodic.kcut;
ewaldParams.boundary = cfg.periodic.boundary;

targetMask = logical(polsys.site_is_polarizable(:));
sourceMask = abs(polsys.site_charge(:)) > 0;

%% ------------------------------------------------------------------------
% Build the two external fields

% nonperiodic external field
params_nonper = struct();
params_nonper.use_thole = cfg.use_thole;
params_nonper.field = struct();
params_nonper.field.mode = 'nonperiodic';
params_nonper.field.exclude_self = true;
params_nonper.field.use_thole_damping = cfg.use_thole;
params_nonper.field.target_mask = targetMask;
params_nonper.field.source_mask = sourceMask;

fprintf('\nBuilding nonperiodic Eext...\n');
Eext_nonper = calc.compute_external_field(polsys, params_nonper);

% periodic external field
params_per = struct();
params_per.use_thole = cfg.use_thole;
params_per.field = struct();
params_per.field.mode = 'periodic';
params_per.field.exclude_self = true;
params_per.field.use_thole_damping = cfg.use_thole;
params_per.field.target_mask = targetMask;
params_per.field.source_mask = sourceMask;
params_per.field.ewald = ewaldParams;
params_per.field.k_block_size = cfg.k_block_size;

fprintf('Building periodic Eext...\n');
Eext_per = calc.compute_external_field(polsys, params_per);

%% ------------------------------------------------------------------------
% Shared solver params

scfParams = struct();
scfParams.use_thole = cfg.use_thole;
scfParams.softening = cfg.softening;
scfParams.verbose = true;
scfParams.tol = cfg.sor.tol;
scfParams.maxIter = cfg.sor.maxIter;
scfParams.omega = cfg.sor.omega;
scfParams.printEvery = cfg.sor.printEvery;
scfParams.residualEvery = cfg.sor.residualEvery;
scfParams.stopMetric = cfg.sor.stopMetric;
scfParams.kspace_mode = cfg.kspace_mode;
scfParams.k_block_size = cfg.k_block_size;
scfParams.kspace_memory_limit_gb = cfg.kspace_memory_limit_gb;

%% ------------------------------------------------------------------------
% Solve case 1: nonperiodic Eext + periodic solver

fprintf('\nPreparing/solving case 1: nonperiodic Eext + periodic SOR\n');
problem_nonper = thole.prepare_scf_problem(polsys, Eext_nonper, scfParams);

cacheParams = struct();
cacheParams.use_thole = cfg.use_thole;
cacheParams.verbose = true;

realCache_nonper = geom.build_periodic_realspace_cache( ...
    polsys, problem_nonper, ewaldParams, cacheParams);

rowCache_nonper = geom.build_periodic_realspace_row_cache( ...
    polsys, problem_nonper, realCache_nonper);

kCache_nonper = geom.build_periodic_kspace_cache( ...
    polsys, problem_nonper, ewaldParams, scfParams);

scfParams_nonper = scfParams;
scfParams_nonper.realspace_cache = realCache_nonper;
scfParams_nonper.realspace_row_cache = rowCache_nonper;
scfParams_nonper.kspace_cache = kCache_nonper;

[mu_nonper, scf_nonper] = thole.solve_scf_iterative_periodic_sor( ...
    polsys, Eext_nonper, ewaldParams, scfParams_nonper);

fprintf('Case 1 done: nIter = %d, relres = %.3e\n', ...
    scf_nonper.nIter, scf_nonper.relres);

%% ------------------------------------------------------------------------
% Solve case 2: periodic Eext + periodic solver

fprintf('\nPreparing/solving case 2: periodic Eext + periodic SOR\n');
problem_per = thole.prepare_scf_problem(polsys, Eext_per, scfParams);

realCache_per = geom.build_periodic_realspace_cache( ...
    polsys, problem_per, ewaldParams, cacheParams);

rowCache_per = geom.build_periodic_realspace_row_cache( ...
    polsys, problem_per, realCache_per);

kCache_per = geom.build_periodic_kspace_cache( ...
    polsys, problem_per, ewaldParams, scfParams);

scfParams_per = scfParams;
scfParams_per.realspace_cache = realCache_per;
scfParams_per.realspace_row_cache = rowCache_per;
scfParams_per.kspace_cache = kCache_per;

[mu_per, scf_per] = thole.solve_scf_iterative_periodic_sor( ...
    polsys, Eext_per, ewaldParams, scfParams_per);

fprintf('Case 2 done: nIter = %d, relres = %.3e\n', ...
    scf_per.nIter, scf_per.relres);

%% ------------------------------------------------------------------------
% Optional diagnostics

local_print_dipole_summary('nonperiodic Eext', polsys, mu_nonper);
local_print_dipole_summary('periodic Eext',    polsys, mu_per);

%% ------------------------------------------------------------------------
% Plot side-by-side

figure('Color', 'w', 'Position', [100 100 1500 700]);

subplot(1,2,1);
local_plot_dipoles(sys, mu_nonper, refID, nbrID, ...
    sprintf('nonperiodic Eext | E = %.3f meV', ...
    1000 * 27.211386245988 * (-0.5 * sum(sum(mu_nonper(targetMask,:) .* Eext_nonper(targetMask,:))))), cfg);

subplot(1,2,2);
local_plot_dipoles(sys, mu_per, refID, nbrID, ...
    sprintf('periodic Eext | E = %.3f meV', ...
    1000 * 27.211386245988 * (-0.5 * sum(sum(mu_per(targetMask,:) .* Eext_per(targetMask,:))))), cfg);

%% ========================================================================
% local helpers
%% ========================================================================

function local_plot_dipoles(sys, mu, refID, nbrID, plotTitle, cfg)

    viz.plot_supercell_selection(sys, ...
        'ReferenceMolID', refID, ...
        'NeighborMolID', nbrID, ...
        'Axes', gca, ...
        'Title', plotTitle, ...
        'DrawBox', true, ...
        'ShowCOM', true, ...
        'ShowLabels', true, ...
        'ShowPairLine', true);

    hold on;

    pos = sys.site_pos;
    muMag = sqrt(sum(mu.^2, 2));

    if cfg.plot.onlyPolarizable
        siteMask = logical(sys.site_is_polarizable(:));
    else
        siteMask = true(size(pos,1),1);
    end

    arrowMask = siteMask & (muMag > cfg.plot.dipoleThreshold);

    if nnz(arrowMask) > cfg.plot.maxArrows
        idx = find(arrowMask);
        [~, order] = maxk(muMag(idx), cfg.plot.maxArrows);
        keep = false(size(arrowMask));
        keep(idx(order)) = true;
        arrowMask = keep;
    end

    p = pos(arrowMask, :);
    u = mu(arrowMask, :);

    quiver3(p(:,1), p(:,2), p(:,3), ...
            cfg.plot.arrowScale * u(:,1), ...
            cfg.plot.arrowScale * u(:,2), ...
            cfg.plot.arrowScale * u(:,3), ...
            0, 'LineWidth', 1.0, 'Color', [0 0 0]);

    hold off;
end

function local_print_dipole_summary(label, polsys, mu)
    q = polsys.site_charge(:);
    r = polsys.site_pos;

    p_src = sum(q .* r, 1);
    p_ind = sum(mu, 1);
    p_tot = p_src + p_ind;

    convD = 2.541746;

    fprintf('\n=== %s ===\n', label);
    fprintf('Bare source |p| = %.6f e*bohr = %.6f D\n', norm(p_src), norm(p_src)*convD);
    fprintf('Induced    |p| = %.6f e*bohr = %.6f D\n', norm(p_ind), norm(p_ind)*convD);
    fprintf('Total      |p| = %.6f e*bohr = %.6f D\n', norm(p_tot), norm(p_tot)*convD);
end

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('Missing direct lattice on system.');
    end
end