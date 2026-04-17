%% iterative_cache_test
clear; clc; close all;

fprintf('=== iterative cache test ===\n');

rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
filename   = fullfile(rootFolder, 'a_0.0_CONTCAR.vasp');

%% ------------------------------------------------------------------------
% User controls

supercellSize = [2 5 1];
bondScale     = 1.20;

% Optional nonperiodic cutoff in bohr. Set [] for exact all-pairs.
rcut_bohr = 15.0;   % [] for exact all-pairs

pairCharges = [+1 -1];

jacobiMixing  = 0.60;
jacobiMaxIter = 200;

sorOmega   = 0.90;
sorMaxIter = 200;

softening = 0.0;

% Comparison tolerances
fieldTol = 1e-12;
muTol    = 1e-8;
eTolHa   = 1e-10;

% Timing controls
nWarmup = 1;
nRepeat = 3;

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
% Build working supercell

opts = struct();
opts.supercell_size = supercellSize;
opts.bondScale      = bondScale;
opts.verbose        = true;
opts.removeMolIDs   = [];
opts.activeMolIDs   = [];

sys = builder.make_crystal_system(crystal, model, opts);

fprintf('\nBuilt working system:\n');
fprintf('  n_unit_sites = %d\n', sys.n_unit_sites);
fprintf('  n_cells      = %d\n', sys.n_cells);
fprintf('  n_sites      = %d\n', sys.n_sites);

%% ------------------------------------------------------------------------
% Automatically choose charged pair

[refID, refSummary] = builder.choose_center_reference_molecule(sys, ...
    'RequireComplete', true, ...
    'Verbose', true);

desc = builder.complete_molecule_descriptors_relative_to_reference(sys, refID, ...
    'Verbose', true);

sameStack = builder.select_same_stack_neighbor(desc, 1, ...
    'Direction', 'either', ...
    'Verbose', true);

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
    'TotalCharges', pairCharges, ...
    'SetActive', true, ...
    'DisablePolarizabilityOnCharged', true, ...
    'ZeroExistingCharges', true, ...
    'Verbose', true);

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
% Canonical AU system

params = struct();
params.ewald = struct();
params.ewald.mode = 'nonperiodic';

polsys = builder.extract_polarization_system(sys, params);
io.assert_atomic_units(polsys);

fprintf('\nCanonical polarization system:\n');
fprintf('  n_sites             = %d\n', polsys.n_sites);
fprintf('  n_polarizable_sites = %d\n', nnz(polsys.site_is_polarizable));

%% ------------------------------------------------------------------------
% External field from charges

fieldParams = struct();
fieldParams.field = struct();
fieldParams.field.include_external_charges = true;
fieldParams.field.exclude_self = true;
fieldParams.field.softening = softening;
if ~isempty(rcut_bohr)
    fieldParams.field.rcut = rcut_bohr;
end

% Build charge-field cache once and thread it into compute_external_field.
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

Eext = calc.compute_external_field(polsys, fieldParams);

if isempty(Eext) || ~isequal(size(Eext), [polsys.n_sites, 3])
    error('External field was not built correctly.');
end
if any(~isfinite(Eext(:)))
    error('External field contains non-finite values.');
end

%% ------------------------------------------------------------------------
% Solver settings: matrix-free Jacobi / damped fixed point

scf = struct();
scf.softening      = softening;
scf.tol            = 1e-8;
scf.maxIter        = jacobiMaxIter;
scf.mixing         = jacobiMixing;
scf.verbose        = true;
scf.printEvery     = 5;
scf.residualEvery  = 10;
scf.stopMetric     = 'max_dmu';

if ~isempty(rcut_bohr)
    scf.rcut = rcut_bohr;
end

problem = thole.prepare_scf_problem(polsys, Eext, scf);

%% ------------------------------------------------------------------------
% Build cache explicitly for A/B checks

cacheOpts = struct();
cacheOpts.site_mask = problem.polMask;
if isfield(scf, 'rcut') && ~isempty(scf.rcut)
    cacheOpts.rcut = scf.rcut;
else
    cacheOpts.rcut = inf;
end

geomCache = geom.build_nonperiodic_pair_cache(polsys, cacheOpts);

fprintf('\nCache diagnostics:\n');
fprintf('  cached site count   = %d\n', numel(geomCache.site_idx));
fprintf('  cached pair count   = %d\n', geomCache.nPairs);
if isfinite(geomCache.rcut)
    fprintf('  cache rcut          = %.6f bohr\n', geomCache.rcut);
else
    fprintf('  cache rcut          = exact (no cutoff)\n');
end

%% ------------------------------------------------------------------------
% Build dipoleParams for direct field comparison

dipoleParams_uncached = struct();
dipoleParams_uncached.exclude_self = true;
dipoleParams_uncached.softening    = softening;
dipoleParams_uncached.target_mask  = problem.polMask;
dipoleParams_uncached.source_mask  = problem.polMask;
if isfield(scf, 'rcut') && ~isempty(scf.rcut)
    dipoleParams_uncached.rcut = scf.rcut;
end

dipoleParams_cached = dipoleParams_uncached;
dipoleParams_cached.geom_cache = geomCache;

%% ------------------------------------------------------------------------
% Compare one induced-field application on a random dipole guess

rng(1);
mu_probe = zeros(polsys.n_sites, 3);
mu_probe(problem.polMask, :) = 1e-3 * randn(nnz(problem.polMask), 3);

Edip_uncached = thole.induced_field_from_dipoles_thole(polsys, mu_probe, dipoleParams_uncached);
Edip_cached   = thole.induced_field_from_dipoles_thole(polsys, mu_probe, dipoleParams_cached);

fieldDiff = max(sqrt(sum((Edip_uncached - Edip_cached).^2, 2)));

fprintf('\nSingle-apply comparison:\n');
fprintf('  max |Edip_uncached - Edip_cached| = %.16e\n', fieldDiff);

if ~isfinite(fieldDiff) || fieldDiff > fieldTol
    error('Cached and uncached induced fields disagree: %.3e', fieldDiff);
end

%% ------------------------------------------------------------------------
% Time repeated field applications

for k = 1:nWarmup
    thole.induced_field_from_dipoles_thole(polsys, mu_probe, dipoleParams_uncached);
    thole.induced_field_from_dipoles_thole(polsys, mu_probe, dipoleParams_cached);
end

t_uncached = zeros(nRepeat, 1);
t_cached   = zeros(nRepeat, 1);

for k = 1:nRepeat
    tic;
    thole.induced_field_from_dipoles_thole(polsys, mu_probe, dipoleParams_uncached);
    t_uncached(k) = toc;

    tic;
    thole.induced_field_from_dipoles_thole(polsys, mu_probe, dipoleParams_cached);
    t_cached(k) = toc;
end

fprintf('\nField-apply timing:\n');
fprintf('  uncached mean = %.6f s\n', mean(t_uncached));
fprintf('  cached   mean = %.6f s\n', mean(t_cached));
fprintf('  speedup       = %.3fx\n', mean(t_uncached) / mean(t_cached));

%% ------------------------------------------------------------------------
% Build explicit operator for energy / relres checks

scfAssembly = scf;
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
% Matrix-free Jacobi / damped fixed-point solve

fprintf('\n--- matrix-free iterative solve ---\n');
tic;
[mu_iter, scf_iter] = thole.solve_scf_iterative(polsys, Eext, scf);
time_iter = toc;

problemIter = thole.prepare_scf_problem(polsys, Eext, scf);
E_iter = calc.compute_total_energy_active_space(polsys, problemIter, mu_iter, Eext, Tpol);
relres_iter = thole.compute_active_space_relres(problemIter, Tpol, mu_iter);

fprintf('iter: relres = %.3e | E = %+0.12f Ha | iters = %d | time = %.3f s | converged = %d\n', ...
    relres_iter, E_iter.total, scf_iter.nIter, time_iter, scf_iter.converged);

%% ------------------------------------------------------------------------
% Matrix-free SOR solve

scf_sor = struct();
scf_sor.softening      = softening;
scf_sor.tol            = 1e-8;
scf_sor.maxIter        = sorMaxIter;
scf_sor.omega          = sorOmega;
scf_sor.verbose        = true;
scf_sor.printEvery     = 5;
scf_sor.residualEvery  = 10;
scf_sor.stopMetric     = 'max_dmu';

if ~isempty(rcut_bohr)
    scf_sor.rcut = rcut_bohr;
end

fprintf('\n--- matrix-free SOR solve ---\n');
tic;
[mu_sor_mf, scf_sor_mf] = thole.solve_scf_iterative_sor(polsys, Eext, scf_sor);
time_sor_mf = toc;

problemSOR = thole.prepare_scf_problem(polsys, Eext, scf_sor);
E_sor_mf = calc.compute_total_energy_active_space(polsys, problemSOR, mu_sor_mf, Eext, Tpol);
relres_sor_mf = thole.compute_active_space_relres(problemSOR, Tpol, mu_sor_mf);

fprintf('iter-SOR: relres = %.3e | E = %+0.12f Ha | iters = %d | time = %.3f s | converged = %d\n', ...
    relres_sor_mf, E_sor_mf.total, scf_sor_mf.nIter, time_sor_mf, scf_sor_mf.converged);

%% ------------------------------------------------------------------------
% Comparisons

fprintf('\n--- comparisons ---\n');

mu_diff_iter = max(sqrt(sum((mu_direct - mu_iter).^2, 2)));
e_diff_iter  = abs(E_direct.total - E_iter.total);

mu_diff_sor = max(sqrt(sum((mu_direct - mu_sor_mf).^2, 2)));
e_diff_sor  = abs(E_direct.total - E_sor_mf.total);

fprintf('Jacobi-like matrix-free:\n');
fprintf('  max |mu_direct - mu_iter| = %.16e\n', mu_diff_iter);
fprintf('  |E_direct - E_iter|       = %.16e Ha\n', e_diff_iter);

fprintf('Matrix-free SOR:\n');
fprintf('  max |mu_direct - mu_sor|  = %.16e\n', mu_diff_sor);
fprintf('  |E_direct - E_sor|        = %.16e Ha\n', e_diff_sor);

%% ------------------------------------------------------------------------
% Assertions

if ~isfinite(direct.relres) || direct.relres > 1e-10
    error('Direct solver relative residual too large: %.3e', direct.relres);
end

if ~scf_iter.converged
    error('Matrix-free iterative solver did not converge.');
end

if ~isfinite(relres_iter) || relres_iter > 1e-6
    error('Matrix-free iterative relative residual too large: %.3e', relres_iter);
end

if mu_diff_iter > muTol
    error('Direct and matrix-free iterative dipoles disagree: %.3e', mu_diff_iter);
end

if e_diff_iter > eTolHa
    error('Direct and matrix-free iterative energies disagree: %.3e Ha', e_diff_iter);
end

if ~scf_sor_mf.converged
    error('Matrix-free SOR solver did not converge.');
end

if ~isfinite(relres_sor_mf) || relres_sor_mf > 1e-6
    error('Matrix-free SOR relative residual too large: %.3e', relres_sor_mf);
end

if mu_diff_sor > muTol
    error('Direct and matrix-free SOR dipoles disagree: %.3e', mu_diff_sor);
end

if e_diff_sor > eTolHa
    error('Direct and matrix-free SOR energies disagree: %.3e Ha', e_diff_sor);
end

fprintf('\nPASS: cached matrix-free Jacobi and SOR solvers match direct reference.\n');