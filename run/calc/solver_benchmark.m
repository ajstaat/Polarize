%% solver_benchmark_nonperiodic
% Benchmark and regression test for nonperiodic polarization solvers.
%
% This script:
%   1) builds a charged dimer-in-crystal test case
%   2) validates cached vs uncached dipole-field application
%   3) benchmarks cached dipole-field application speed
%   4) assembles the nonperiodic active-space operator
%   5) compares:
%        - direct solve
%        - matrix-free iterative (Jacobi-like / damped fixed-point)
%        - matrix-free iterative SOR
%
% Intended use:
%   - regression check after solver/kernel changes
%   - benchmark for nonperiodic matrix-free performance
%
% Notes:
%   - tolerances below are tuned for the current fast-stop benchmark mode
%   - matrix-free iterative methods are compared against the direct reference

clear; clc; close all;

fprintf('=== nonperiodic solver benchmark ===\n');

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

cfg.jacobi = struct();
cfg.jacobi.mixing        = 0.60;
cfg.jacobi.maxIter       = 200;
cfg.jacobi.tol           = 1e-8;
cfg.jacobi.printEvery    = 5;
cfg.jacobi.residualEvery = 10;
cfg.jacobi.stopMetric    = 'max_dmu';

cfg.sor = struct();
cfg.sor.omega         = 0.90;
cfg.sor.maxIter       = 200;
cfg.sor.tol           = 1e-8;
cfg.sor.printEvery    = 5;
cfg.sor.residualEvery = 10;
cfg.sor.stopMetric    = 'max_dmu';

cfg.verbose = true;

% Timing controls for repeated field-apply benchmark
cfg.nWarmup = 1;
cfg.nRepeat = 3;

% Regression tolerances
tol = struct();
tol.field = 1e-12;
tol.mu    = 1e-8;
tol.eHa   = 1e-10;
tol.relres_iter = 1e-6;
tol.relres_sor  = 1e-6;
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
% External field from charges (with cache)

fieldParams = struct();
fieldParams.field = struct();
fieldParams.field.include_external_charges = true;
fieldParams.field.exclude_self = true;
fieldParams.field.softening = cfg.softening;
if ~isempty(cfg.rcut_bohr)
    fieldParams.field.rcut = cfg.rcut_bohr;
end

fieldCacheOpts = struct();
fieldCacheOpts.site_mask = true(polsys.n_sites, 1);

if isfield(fieldParams.field, 'target_mask') && ~isempty(fieldParams.field.target_mask) && ...
   isfield(fieldParams.field, 'source_mask') && ~isempty(fieldParams.field.source_mask)
    fieldCacheOpts.site_mask = logical(fieldParams.field.target_mask(:)) | ...
                               logical(fieldParams.field.source_mask(:));
elseif isfield(fieldParams.field, 'target_mask') && ~isempty(fieldParams.field.target_mask)
    fieldCacheOpts.site_mask = logical(fieldParams.field.target_mask(:));
elseif isfield(fieldParams.field, 'source_mask') && ~isempty(fieldParams.field.source_mask)
    fieldCacheOpts.site_mask = logical(fieldParams.field.source_mask(:));
end

if isfield(fieldParams.field, 'rcut') && ~isempty(fieldParams.field.rcut)
    fieldCacheOpts.rcut = fieldParams.field.rcut;
else
    fieldCacheOpts.rcut = inf;
end

fieldParams.field.geom_cache = geom.build_nonperiodic_pair_cache(polsys, fieldCacheOpts);

tField = tic;
Eext = calc.compute_external_field(polsys, fieldParams);
time_external_field = toc(tField);

if isempty(Eext) || ~isequal(size(Eext), [polsys.n_sites, 3])
    error('External field was not built correctly.');
end
if any(~isfinite(Eext(:)))
    error('External field contains non-finite values.');
end

fprintf('\nExternal field:\n');
fprintf('  build time = %.6f s\n', time_external_field);

%% ------------------------------------------------------------------------
% Baseline SCF config for active-space bookkeeping and dipole cache checks

scfBase = struct();
scfBase.softening      = cfg.softening;
scfBase.tol            = cfg.jacobi.tol;
scfBase.maxIter        = cfg.jacobi.maxIter;
scfBase.mixing         = cfg.jacobi.mixing;
scfBase.verbose        = cfg.verbose;
scfBase.printEvery     = cfg.jacobi.printEvery;
scfBase.residualEvery  = cfg.jacobi.residualEvery;
scfBase.stopMetric     = cfg.jacobi.stopMetric;
if ~isempty(cfg.rcut_bohr)
    scfBase.rcut = cfg.rcut_bohr;
end

problem = thole.prepare_scf_problem(polsys, Eext, scfBase);

%% ------------------------------------------------------------------------
% Dipole pair cache diagnostics

cacheOpts = struct();
cacheOpts.site_mask = problem.polMask;
if isfield(scfBase, 'rcut') && ~isempty(scfBase.rcut)
    cacheOpts.rcut = scfBase.rcut;
else
    cacheOpts.rcut = inf;
end

geomCache = geom.build_nonperiodic_pair_cache(polsys, cacheOpts);

fprintf('\nDipole cache diagnostics:\n');
fprintf('  cached site count   = %d\n', numel(geomCache.site_idx));
fprintf('  cached pair count   = %d\n', geomCache.nPairs);
if isfinite(geomCache.rcut)
    fprintf('  cache rcut          = %.6f bohr\n', geomCache.rcut);
else
    fprintf('  cache rcut          = exact (no cutoff)\n');
end

%% ------------------------------------------------------------------------
% Validate cached vs uncached dipole-field application

dipoleParams_uncached = struct();
dipoleParams_uncached.exclude_self = true;
dipoleParams_uncached.softening    = cfg.softening;
dipoleParams_uncached.target_mask  = problem.polMask;
dipoleParams_uncached.source_mask  = problem.polMask;
if isfield(scfBase, 'rcut') && ~isempty(scfBase.rcut)
    dipoleParams_uncached.rcut = scfBase.rcut;
end

dipoleParams_cached = dipoleParams_uncached;
dipoleParams_cached.geom_cache = geomCache;

rng(1);
mu_probe = zeros(polsys.n_sites, 3);
mu_probe(problem.polMask, :) = 1e-3 * randn(nnz(problem.polMask), 3);

Edip_uncached = thole.induced_field_from_dipoles_thole(polsys, mu_probe, dipoleParams_uncached);
Edip_cached   = thole.induced_field_from_dipoles_thole(polsys, mu_probe, dipoleParams_cached);

fieldDiff = max(sqrt(sum((Edip_uncached - Edip_cached).^2, 2)));

fprintf('\nSingle-apply validation:\n');
fprintf('  max |Edip_uncached - Edip_cached| = %.16e\n', fieldDiff);

if ~isfinite(fieldDiff) || fieldDiff > tol.field
    error('Cached and uncached induced fields disagree: %.3e', fieldDiff);
end

%% ------------------------------------------------------------------------
% Benchmark repeated dipole-field applications

for k = 1:cfg.nWarmup
    thole.induced_field_from_dipoles_thole(polsys, mu_probe, dipoleParams_uncached);
    thole.induced_field_from_dipoles_thole(polsys, mu_probe, dipoleParams_cached);
end

t_uncached = zeros(cfg.nRepeat, 1);
t_cached   = zeros(cfg.nRepeat, 1);

for k = 1:cfg.nRepeat
    tic;
    thole.induced_field_from_dipoles_thole(polsys, mu_probe, dipoleParams_uncached);
    t_uncached(k) = toc;

    tic;
    thole.induced_field_from_dipoles_thole(polsys, mu_probe, dipoleParams_cached);
    t_cached(k) = toc;
end

fprintf('\nDipole field-apply timing:\n');
fprintf('  uncached mean = %.6f s\n', mean(t_uncached));
fprintf('  cached   mean = %.6f s\n', mean(t_cached));
fprintf('  speedup       = %.3fx\n', mean(t_uncached) / mean(t_cached));

%% ------------------------------------------------------------------------
% Assemble active-space operator for direct reference and residual checks

scfAssembly = scfBase;
scfAssembly.verbose = true;
scfAssembly.printEvery = 250;
scfAssembly.residualEvery = 250;

[Tpol, opinfo] = ewald.assemble_nonperiodic_interaction_matrix(polsys, problem, struct(), scfAssembly);

fprintf('\nOperator diagnostics:\n');
fprintf('  active-space operator size = %d x %d\n', size(Tpol,1), size(Tpol,2));
fprintf('  nPolSites                  = %d\n', problem.nPolSites);
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
% Matrix-free iterative solve (Jacobi-like / damped fixed point)

scf_iter = struct();
scf_iter.softening      = cfg.softening;
scf_iter.tol            = cfg.jacobi.tol;
scf_iter.maxIter        = cfg.jacobi.maxIter;
scf_iter.mixing         = cfg.jacobi.mixing;
scf_iter.verbose        = cfg.verbose;
scf_iter.printEvery     = cfg.jacobi.printEvery;
scf_iter.residualEvery  = cfg.jacobi.residualEvery;
scf_iter.stopMetric     = cfg.jacobi.stopMetric;
if ~isempty(cfg.rcut_bohr)
    scf_iter.rcut = cfg.rcut_bohr;
end

fprintf('\n--- matrix-free iterative solve ---\n');
tic;
[mu_iter, scf_iter_out] = thole.solve_scf_iterative(polsys, Eext, scf_iter);
time_iter = toc;

problemIter = thole.prepare_scf_problem(polsys, Eext, scf_iter);
E_iter = calc.compute_total_energy_active_space(polsys, problemIter, mu_iter, Eext, Tpol);
relres_iter = thole.compute_active_space_relres(problemIter, Tpol, mu_iter);

fprintf('iter: relres = %.3e | E = %+0.12f Ha | iters = %d | time = %.3f s | converged = %d\n', ...
    relres_iter, E_iter.total, scf_iter_out.nIter, time_iter, scf_iter_out.converged);

%% ------------------------------------------------------------------------
% Matrix-free SOR solve

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

mu_diff_iter = max(sqrt(sum((mu_direct - mu_iter).^2, 2)));
e_diff_iter  = abs(E_direct.total - E_iter.total);

mu_diff_sor = max(sqrt(sum((mu_direct - mu_sor).^2, 2)));
e_diff_sor  = abs(E_direct.total - E_sor.total);

fprintf('Jacobi-like matrix-free:\n');
fprintf('  max |mu_direct - mu_iter| = %.16e\n', mu_diff_iter);
fprintf('  |E_direct - E_iter|       = %.16e Ha\n', e_diff_iter);

fprintf('Matrix-free SOR:\n');
fprintf('  max |mu_direct - mu_sor|  = %.16e\n', mu_diff_sor);
fprintf('  |E_direct - E_sor|        = %.16e Ha\n', e_diff_sor);

fprintf('\n--- timing summary ---\n');
fprintf('  external field        = %.6f s\n', time_external_field);
fprintf('  direct solve          = %.6f s\n', time_direct);
fprintf('  matrix-free iterative = %.6f s\n', time_iter);
fprintf('  matrix-free SOR       = %.6f s\n', time_sor);

%% ------------------------------------------------------------------------
% Assertions

if ~isfinite(direct.relres) || direct.relres > tol.relres_direct
    error('Direct solver relative residual too large: %.3e', direct.relres);
end

if ~scf_iter_out.converged
    error('Matrix-free iterative solver did not converge.');
end
if ~isfinite(relres_iter) || relres_iter > tol.relres_iter
    error('Matrix-free iterative relative residual too large: %.3e', relres_iter);
end
if mu_diff_iter > tol.mu
    error('Direct and matrix-free iterative dipoles disagree: %.3e', mu_diff_iter);
end
if e_diff_iter > tol.eHa
    error('Direct and matrix-free iterative energies disagree: %.3e Ha', e_diff_iter);
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

fprintf('\nPASS: nonperiodic solver benchmark completed successfully.\n');