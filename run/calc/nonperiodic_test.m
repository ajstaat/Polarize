%% nonperiodic_real_crystal_example
% Example: nonperiodic polarization calculation on a real charged-pair crystal.

clear; clc; close all;

fprintf('=== nonperiodic polarization energy (real crystal) ===\n');

%% ------------------------------------------------------------------------
% User controls

cfg = struct();

cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename   = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

cfg.supercellSize = [2 5 1];
cfg.bondScale     = 1.20;

cfg.pairCharges = [+1 -1];

cfg.use_thole = true;
cfg.verbose = true;
cfg.softening = 0.0;

cfg.nonperiodic = struct();
cfg.nonperiodic.softening = 0.0;
cfg.nonperiodic.use_thole = cfg.use_thole;
cfg.nonperiodic.rcut = 18.0;

%% ------------------------------------------------------------------------
% Import crystal template

fprintf('\n[setup] importing crystal template...\n');
crystal = io.import_contcar_as_crystal(cfg.filename, ...
    'BondScale', cfg.bondScale, ...
    'WrapFractional', true, ...
    'SortMolecules', false);

fprintf('Imported crystal template:\n');
fprintf('  nSites     = %d\n', size(crystal.cart_coords, 1));
fprintf('  nBaseMols  = %d\n', numel(unique(crystal.base_mol_id)));

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

fprintf('\n[setup] building crystal system...\n');
buildOpts = struct();
buildOpts.supercell_size = cfg.supercellSize;
buildOpts.bondScale      = cfg.bondScale;
buildOpts.verbose        = cfg.verbose;
buildOpts.removeMolIDs   = [];
buildOpts.activeMolIDs   = [];

sys0 = builder.make_crystal_system(crystal, model, buildOpts);

fprintf('\nBuilt working system:\n');
fprintf('  n_unit_sites = %d\n', sys0.n_unit_sites);
fprintf('  n_cells      = %d\n', sys0.n_cells);
fprintf('  n_sites      = %d\n', sys0.n_sites);

%% ------------------------------------------------------------------------
% Automatically choose charged pair

fprintf('\n[setup] choosing charged-pair molecules...\n');
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

fprintf('\nChosen pair (automatic):\n');
fprintf('  ref ID = %d\n', refID);
fprintf('  nbr ID = %d\n', nbrID);

%% ------------------------------------------------------------------------
% Apply charges and disable polarizability on charged molecules

fprintf('\n[setup] applying charges and extracting polarization system...\n');
sys = builder.apply_molecule_charges(sys0, [refID nbrID], ...
    'Mode', 'uniform', ...
    'TotalCharges', cfg.pairCharges, ...
    'SetActive', true, ...
    'DisablePolarizabilityOnCharged', true, ...
    'ZeroExistingCharges', true, ...
    'Verbose', cfg.verbose);

params_extract = struct();
params_extract.ewald = struct();
params_extract.ewald.mode = 'periodic';  % extraction-path consistency

polsys0 = builder.extract_polarization_system(sys, params_extract);
io.assert_atomic_units(polsys0);

H0 = local_get_direct_lattice(polsys0);
polsys = polsys0;
polsys.super_lattice = H0;

center0 = 0.5 * sum(H0, 2);
center_per = 0.5 * sum(polsys.super_lattice, 2);
shift = (center_per - center0).';
polsys.site_pos = polsys0.site_pos + shift;

fprintf('\nCanonical polarization system:\n');
fprintf('  n_sites             = %d\n', polsys.n_sites);
fprintf('  n_polarizable_sites = %d\n', nnz(polsys.site_is_polarizable));

targetMask = logical(polsys.site_is_polarizable(:));
sourceMask = abs(polsys.site_charge(:)) > 0;

fprintf('\nMask summary:\n');
fprintf('  n target sites (polarizable) = %d\n', nnz(targetMask));
fprintf('  n source sites (charged)     = %d\n', nnz(sourceMask));
fprintf('  total source charge          = %+0.12f\n', sum(polsys.site_charge(sourceMask)));
fprintf('  max |source charge|          = %.12e\n', max(abs(polsys.site_charge(sourceMask))));

%% ------------------------------------------------------------------------
% Nonperiodic SOR controls

sorParams = struct();
sorParams.use_thole = cfg.use_thole;
sorParams.softening = cfg.softening;
sorParams.verbose = true;
sorParams.printEvery = 10;
sorParams.residualEvery = 10;
sorParams.tol = 1e-6;
sorParams.maxIter = 500;
sorParams.omega = 0.97;
sorParams.stopMetric = 'max_dmu';
sorParams.rcut = cfg.nonperiodic.rcut;
sorParams.checkResidualAgainstLegacy = false;

%% ------------------------------------------------------------------------
% NONPERIODIC run
fprintf('\n============================================================\n');
fprintf('NONPERIODIC\n');
fprintf('============================================================\n');

params_nonper = struct();
params_nonper.use_thole = cfg.use_thole;
params_nonper.field = struct();
params_nonper.field.mode = 'nonperiodic';
params_nonper.field.exclude_self = true;
params_nonper.field.use_thole_damping = cfg.use_thole;
params_nonper.field.target_mask = targetMask;
params_nonper.field.source_mask = sourceMask;
% Intentionally no field rcut here: preserve old uncut external field.
% Intentionally no field geom_cache here: use optimized direct evaluator.

fprintf('[nonperiodic] building external field...\n');
tField = tic;
Eext = calc.compute_external_field(polsys, params_nonper);
timeField = toc(tField);
fprintf('[nonperiodic] external field done in %.6f s | ||Eext||_F = %.16e\n', ...
    timeField, norm(Eext, 'fro'));

fprintf('[nonperiodic] preparing SCF problem...\n');
tProblem = tic;
problem = thole.prepare_scf_problem(polsys, Eext, sorParams);
timeProblem = toc(tProblem);
fprintf('[nonperiodic] problem prep done in %.6f s | nPolSites = %d\n', ...
    timeProblem, problem.nPolSites);

fprintf('[nonperiodic] building solver row cache...\n');
tRowCache = tic;
rowOpts = struct();
rowOpts.rcut = cfg.nonperiodic.rcut;
rowOpts.profile = true;
rowOpts.use_mex = true;

% Prebuild active-space spatial index once and pass it in.
posAct = polsys.site_pos(problem.activeSites, :);
spatialOpts = struct();
spatialOpts.isPeriodic = false;
spatialOpts.method = 'auto';
spatialOpts.cutoff = cfg.nonperiodic.rcut;
spatial = geom.build_spatial_index(posAct, spatialOpts);

solveRowCache = geom.build_active_row_cache(polsys, problem, rowOpts, spatial);
timeRowCache = toc(tRowCache);
fprintf('[nonperiodic] solver row cache done in %.6f s | nActive = %d\n', ...
    timeRowCache, solveRowCache.nActive);

sorParams.row_cache = solveRowCache;
if isfield(sorParams, 'geom_cache')
    sorParams = rmfield(sorParams, 'geom_cache');
end

fprintf('[nonperiodic] solving matrix-free SOR...\n');
tSolve = tic;
[mu, iterInfo] = thole.solve_scf_iterative_sor(polsys, Eext, sorParams);
timeSolve = toc(tSolve);
fprintf('[nonperiodic] matrix-free SOR done in %.6f s | ||mu||_2 = %.16e\n', ...
    timeSolve, norm(mu, 'fro'));
fprintf('  converged = %d | nIter = %d | relres = %.16e\n', ...
    iterInfo.converged, iterInfo.nIter, iterInfo.relres);

polMask = logical(polsys.site_is_polarizable(:));
E_total = -0.5 * sum(sum(mu(polMask, :) .* Eext(polMask, :)));
fprintf('matrix-free total energy = %+0.12e Ha\n', E_total);

fprintf('\nDone.\n');

%% ========================================================================
% Local helpers
%% ========================================================================

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('Missing direct lattice on system.');
    end
end