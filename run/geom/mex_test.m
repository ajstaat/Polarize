%% test_nonperiodic_pair_cache_direct_vs_sor
% Regression / benchmark test for nonperiodic pair-cache acceleration.
%
% This script:
%   1) builds a charged dimer-in-crystal nonperiodic test case
%   2) builds the nonperiodic pair cache with MEX enabled
%   3) assembles the nonperiodic active-space operator
%   4) compares:
%        - direct solve
%        - matrix-free SOR
%
% Intended use:
%   - validate mex_build_nonperiodic_pair_cache path
%   - benchmark pair-cache build timing
%   - confirm direct and SOR still agree numerically
%
% Notes:
%   - uses a modest supercell so dense direct solve remains feasible
%   - leaves external-field construction cached, as in solver benchmark style

clear; clc; close all;

fprintf('=== nonperiodic pair-cache direct-vs-SOR test ===\n');

%% ------------------------------------------------------------------------
% User controls

cfg = struct();

% Input structure
cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename   = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

% Crystal/system construction
cfg.supercellSize = [2 5 1];
cfg.bondScale     = 1.20;

% Nonperiodic cutoff in bohr; set [] for exact all-pairs
cfg.rcut_bohr = 15.0;

% Charged-pair setup
cfg.pairCharges = [+1 -1];

% Solver controls
cfg.softening = 0.0;
cfg.use_thole = true;
cfg.verbose = true;

cfg.sor = struct();
cfg.sor.omega         = 0.90;
cfg.sor.maxIter       = 200;
cfg.sor.tol           = 1e-8;
cfg.sor.printEvery    = 5;
cfg.sor.residualEvery = 10;
cfg.sor.stopMetric    = 'max_dmu';

% Regression tolerances
tol = struct();
tol.mu            = 1e-8;
tol.eHa           = 1e-10;
tol.relres_sor    = 1e-6;
tol.relres_direct = 1e-10;

%% ------------------------------------------------------------------------
% Import crystal template

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

buildOpts = struct();
buildOpts.supercell_size = cfg.supercellSize;
buildOpts.bondScale      = cfg.bondScale;
buildOpts.verbose        = cfg.verbose;
buildOpts.removeMolIDs   = [];
buildOpts.activeMolIDs   = [];

sys = builder.make_crystal_system(crystal, model, buildOpts);

fprintf('\nBuilt working system:\n');
fprintf('  n_unit_sites = %d\n', sys.n_unit_sites);
fprintf('  n_cells      = %d\n', sys.n_cells);
fprintf('  n_sites      = %d\n', sys.n_sites);

%% ------------------------------------------------------------------------
% Automatically choose charged pair

[refID, refSummary] = builder.choose_center_reference_molecule(sys, ...
    'RequireComplete', true, ...
    'Verbose', cfg.verbose);

desc = builder.complete_molecule_descriptors_relative_to_reference(sys, refID, ...
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

iRef = find(refSummary.candidate_mol_ids == refID, 1, 'first');
if ~isempty(iRef)
    fprintf('  ref distance to center = %.4f\n', refSummary.candidate_distance(iRef));
end

%% ------------------------------------------------------------------------
% Apply charges and disable polarizability on charged molecules

sys = builder.apply_molecule_charges(sys, [refID nbrID], ...
    'Mode', 'uniform', ...
    'TotalCharges', cfg.pairCharges, ...
    'SetActive', true, ...
    'DisablePolarizabilityOnCharged', true, ...
    'ZeroExistingCharges', true, ...
    'Verbose', cfg.verbose);

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

%% ------------------------------------------------------------------------
% Canonical atomic-unit polarization system

params = struct();
params.ewald = struct();
params.ewald.mode = 'nonperiodic';

polsys = builder.extract_polarization_system(sys, params);
io.assert_atomic_units(polsys);

fprintf('\nCanonical polarization system:\n');
fprintf('  n_sites             = %d\n', polsys.n_sites);
fprintf('  n_polarizable_sites = %d\n', nnz(polsys.site_is_polarizable));

%% ------------------------------------------------------------------------
% External field from charges (cached)

fieldParams = struct();
fieldParams.field = struct();
fieldParams.field.include_external_charges = true;
fieldParams.field.exclude_self = true;
fieldParams.field.softening = cfg.softening;
fieldParams.field.use_thole_damping = cfg.use_thole;

if ~isempty(cfg.rcut_bohr)
    fieldParams.field.rcut = cfg.rcut_bohr;
end

fieldCacheOpts = struct();
fieldCacheOpts.site_mask = true(polsys.n_sites, 1);
if isfield(fieldParams.field, 'rcut') && ~isempty(fieldParams.field.rcut)
    fieldCacheOpts.rcut = fieldParams.field.rcut;
else
    fieldCacheOpts.rcut = inf;
end
fieldCacheOpts.use_mex = true;
fieldCacheOpts.profile = true;

tFieldCache = tic;
fieldGeomCache = geom.build_nonperiodic_pair_cache(polsys, fieldCacheOpts);
timeFieldCache = toc(tFieldCache);

fieldParams.field.geom_cache = fieldGeomCache;

tField = tic;
Eext = calc.compute_external_field(polsys, fieldParams);
timeField = toc(tField);

if isempty(Eext) || ~isequal(size(Eext), [polsys.n_sites, 3])
    error('External field was not built correctly.');
end
if any(~isfinite(Eext(:)))
    error('External field contains non-finite values.');
end

fprintf('\nExternal field:\n');
fprintf('  pair-cache build time = %.6f s\n', timeFieldCache);
fprintf('  field build time      = %.6f s\n', timeField);

%% ------------------------------------------------------------------------
% Baseline SCF config and active-space bookkeeping

scfBase = struct();
scfBase.softening      = cfg.softening;
scfBase.tol            = cfg.sor.tol;
scfBase.maxIter        = cfg.sor.maxIter;
scfBase.verbose        = cfg.verbose;
scfBase.printEvery     = cfg.sor.printEvery;
scfBase.residualEvery  = cfg.sor.residualEvery;
scfBase.stopMetric     = cfg.sor.stopMetric;
if ~isempty(cfg.rcut_bohr)
    scfBase.rcut = cfg.rcut_bohr;
end

problem = thole.prepare_scf_problem(polsys, Eext, scfBase);

%% ------------------------------------------------------------------------
% Dipole pair cache with MEX enabled

cacheOpts = struct();
cacheOpts.site_mask = problem.polMask;
if isfield(scfBase, 'rcut') && ~isempty(scfBase.rcut)
    cacheOpts.rcut = scfBase.rcut;
else
    cacheOpts.rcut = inf;
end
cacheOpts.use_mex = true;
cacheOpts.profile = true;

% Prebuild spatial index once and pass it in.
posPol = polsys.site_pos(problem.polMask, :);
spatialOpts = struct();
spatialOpts.isPeriodic = false;
if isfinite(cacheOpts.rcut)
    spatialOpts.method = 'auto';
    spatialOpts.cutoff = cacheOpts.rcut;
else
    spatialOpts.method = 'bruteforce';
end

tSpatial = tic;
spatial = geom.build_spatial_index(polsys.site_pos, spatialOpts);
timeSpatial = toc(tSpatial);

tCache = tic;
geomCache = geom.build_nonperiodic_pair_cache(polsys, cacheOpts, spatial);
timeCache = toc(tCache);

fprintf('\nDipole cache diagnostics:\n');
fprintf('  spatial build time    = %.6f s\n', timeSpatial);
fprintf('  pair-cache build time = %.6f s\n', timeCache);
fprintf('  cached site count     = %d\n', numel(geomCache.site_idx));
fprintf('  cached pair count     = %d\n', geomCache.nPairs);
if isfinite(geomCache.rcut)
    fprintf('  cache rcut            = %.6f bohr\n', geomCache.rcut);
else
    fprintf('  cache rcut            = exact (no cutoff)\n');
end

%% ------------------------------------------------------------------------
% Assemble active-space operator for direct reference

scfAssembly = scfBase;
scfAssembly.verbose = true;
scfAssembly.printEvery = 250;
scfAssembly.residualEvery = 250;

tAssemble = tic;
[Tpol, opinfo] = ewald.assemble_nonperiodic_interaction_matrix(polsys, problem, struct(), scfAssembly);
timeAssemble = toc(tAssemble);

fprintf('\nOperator diagnostics:\n');
fprintf('  active-space operator size = %d x %d\n', size(Tpol,1), size(Tpol,2));
fprintf('  nPolSites                  = %d\n', problem.nPolSites);
fprintf('  assembly time              = %.6f s\n', timeAssemble);
fprintf('  opinfo.space               = %s\n', opinfo.space);
if isfield(opinfo, 'rcut') && isfinite(opinfo.rcut)
    nPol = problem.nPolSites;
    nPairsFull = nPol * (nPol - 1) / 2;
    fprintf('  opinfo.rcut                = %.6f bohr\n', opinfo.rcut);
    fprintf('  pair blocks kept           = %d / %d (%.4f%%)\n', ...
        opinfo.nPairBlocksKept, nPairsFull, ...
        100 * opinfo.nPairBlocksKept / nPairsFull);
else
    fprintf('  opinfo.rcut                = exact (no cutoff)\n');
end

%% ------------------------------------------------------------------------
% Direct reference solve

fprintf('\n--- direct solve ---\n');
tic;
[mu_direct, direct] = thole.solve_scf_direct(problem, Tpol);
time_direct = toc;

E_direct = calc.compute_total_energy_active_space(polsys, problem, mu_direct, Eext, Tpol);
relres_direct = thole.compute_active_space_relres(problem, Tpol, mu_direct);

fprintf('direct: relres = %.3e | E = %+0.12f Ha | time = %.3f s\n', ...
    relres_direct, E_direct.total, time_direct);

%% ------------------------------------------------------------------------
% Matrix-free SOR solve with row cache

scf_sor = struct();
scf_sor.softening      = cfg.softening;
scf_sor.tol            = cfg.sor.tol;
scf_sor.maxIter        = cfg.sor.maxIter;
scf_sor.omega          = cfg.sor.omega;
scf_sor.verbose        = cfg.verbose;
scf_sor.printEvery     = cfg.sor.printEvery;
scf_sor.residualEvery  = cfg.sor.residualEvery;
scf_sor.stopMetric     = cfg.sor.stopMetric;
if ~isempty(cfg.rcut_bohr)
    scf_sor.rcut = cfg.rcut_bohr;
end

rowOpts = struct();
rowOpts.rcut = cacheOpts.rcut;
rowOpts.use_mex = true;
rowOpts.profile = true;

posAct = polsys.site_pos(problem.activeSites, :);
rowSpatialOpts = struct();
rowSpatialOpts.isPeriodic = false;
if isfinite(rowOpts.rcut)
    rowSpatialOpts.method = 'auto';
    rowSpatialOpts.cutoff = rowOpts.rcut;
else
    rowSpatialOpts.method = 'bruteforce';
end
rowSpatial = geom.build_spatial_index(posAct, rowSpatialOpts);

tRowCache = tic;
rowCache = geom.build_active_row_cache(polsys, problem, rowOpts, rowSpatial);
timeRowCache = toc(tRowCache);

scf_sor.row_cache = rowCache;

fprintf('\n--- matrix-free SOR solve ---\n');
tic;
[mu_sor, scf_sor_out] = thole.solve_scf_iterative_sor(polsys, Eext, scf_sor);
time_sor = toc;

problemSOR = thole.prepare_scf_problem(polsys, Eext, scf_sor);
E_sor = calc.compute_total_energy_active_space(polsys, problemSOR, mu_sor, Eext, Tpol);
relres_sor = thole.compute_active_space_relres(problemSOR, Tpol, mu_sor);

fprintf('iter-SOR: relres = %.3e | E = %+0.12f Ha | iters = %d | time = %.3f s | converged = %d\n', ...
    relres_sor, E_sor.total, scf_sor_out.nIter, time_sor, scf_sor_out.converged);

%% ------------------------------------------------------------------------
% Comparisons

fprintf('\n--- comparisons ---\n');

mu_diff_sor = max(sqrt(sum((mu_direct - mu_sor).^2, 2)));
e_diff_sor  = abs(E_direct.total - E_sor.total);

fprintf('Matrix-free SOR:\n');
fprintf('  max |mu_direct - mu_sor|  = %.16e\n', mu_diff_sor);
fprintf('  |E_direct - E_sor|        = %.16e Ha\n', e_diff_sor);

fprintf('\n--- timing summary ---\n');
fprintf('  field pair cache      = %.6f s\n', timeFieldCache);
fprintf('  external field        = %.6f s\n', timeField);
fprintf('  dipole pair cache     = %.6f s\n', timeCache);
fprintf('  operator assembly     = %.6f s\n', timeAssemble);
fprintf('  row cache             = %.6f s\n', timeRowCache);
fprintf('  direct solve          = %.6f s\n', time_direct);
fprintf('  matrix-free SOR       = %.6f s\n', time_sor);

%% ------------------------------------------------------------------------
% Assertions

if ~isfinite(direct.relres) || direct.relres > tol.relres_direct
    error('Direct solver relative residual too large: %.3e', direct.relres);
end

if ~scf_sor_out.converged
    error('Matrix-free SOR solver did not converge.');
end
if ~isfinite(relres_sor) || relres_sor > tol.relres_sor
    error('Matrix-free SOR relative residual too large: %.3e', relres_sor);
end
if mu_diff_sor > tol.mu
    error('Direct and matrix-free SOR dipoles disagree: %.3e', mu_diff_sor);
end
if e_diff_sor > tol.eHa
    error('Direct and matrix-free SOR energies disagree: %.3e Ha', e_diff_sor);
end

fprintf('\nPASS: nonperiodic pair-cache direct-vs-SOR test completed successfully.\n');