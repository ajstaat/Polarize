%% periodic_operator_benchmark_active_row_cache_periodic
% Focused benchmark / regression test for periodic polarization pieces
% using the periodic cell-list direct active-row-cache path.
%
% This script:
%   1) builds a charged-pair periodic test case
%   2) builds periodic external field
%   3) prepares the SCF problem
%   4) builds the periodic active row cache (cell-list/MEX path)
%   5) builds the periodic k-space cache
%   6) benchmarks one periodic induced-field application on a probe dipole
%   7) runs periodic SOR using the periodic active row cache
%
% Intended use:
%   - benchmark the periodic cell-list real-space builder
%   - compare directly against the older brute-force periodic row-cache path
%   - keep the rest of the periodic solver path unchanged

clear; clc; close all;

fprintf('=== periodic operator benchmark (cell-list active row cache) ===\n');

%% ------------------------------------------------------------------------
% User controls

cfg = struct();

cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename   = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

cfg.supercellSize = [2 5 1];
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
cfg.sor.maxIter       = 60;
cfg.sor.omega         = 0.45;
cfg.sor.printEvery    = 10;
cfg.sor.residualEvery = 10;
cfg.sor.stopMetric    = 'max_dmu';

cfg.kspace_mode = 'auto';
cfg.k_block_size = 2048;
cfg.kspace_memory_limit_gb = 8;

cfg.nWarmup = 1;
cfg.nRepeat = 3;
cfg.probeDipoleScale = 1e-3;

cfg.profile_row_cache = true;
cfg.use_mex = true;

% Optional: manually override periodic cell-list grid shape
% cfg.grid_shape = [8 8 8];

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
% Choose charged pair automatically

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

fprintf('\nChosen pair:\n');
fprintf('  ref ID = %d\n', refID);
fprintf('  nbr ID = %d\n', nbrID);

%% ------------------------------------------------------------------------
% Apply charges and extract periodic polarization system

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

fprintf('\nCanonical polarization system:\n');
fprintf('  n_sites             = %d\n', polsys.n_sites);
fprintf('  n_polarizable_sites = %d\n', nnz(polsys.site_is_polarizable));

ewaldParams = struct();
ewaldParams.alpha    = cfg.periodic.alpha;
ewaldParams.rcut     = cfg.periodic.rcut;
ewaldParams.kcut     = cfg.periodic.kcut;
ewaldParams.boundary = cfg.periodic.boundary;

targetMask = logical(polsys.site_is_polarizable(:));
sourceMask = abs(polsys.site_charge(:)) > 0;

fprintf('\nMask summary:\n');
fprintf('  n target sites (polarizable) = %d\n', nnz(targetMask));
fprintf('  n source sites (charged)     = %d\n', nnz(sourceMask));
fprintf('  total source charge          = %+0.12f\n', sum(polsys.site_charge(sourceMask)));

%% ------------------------------------------------------------------------
% Periodic external field
%
% Leave this path unchanged for now. The point of this script is the new
% periodic solver-side real-space cache.

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

fprintf('\n[periodic] building Eext...\n');
tEext = tic;
Eext = calc.compute_external_field(polsys, params_per);
timeEext = toc(tEext);

fprintf('[periodic] Eext done in %.6f s | ||Eext||_F = %.16e\n', ...
    timeEext, norm(Eext, 'fro'));

%% ------------------------------------------------------------------------
% Shared solver params / problem

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
scfParams.profile = cfg.profile_row_cache;
scfParams.use_mex = cfg.use_mex;

if isfield(cfg, 'grid_shape')
    scfParams.grid_shape = cfg.grid_shape;
end

fprintf('\n[periodic] preparing SCF problem...\n');
tProblem = tic;
problem = thole.prepare_scf_problem(polsys, Eext, scfParams);
timeProblem = toc(tProblem);

fprintf('[periodic] problem prep done in %.6f s | nPolSites = %d\n', ...
    timeProblem, problem.nPolSites);

%% ------------------------------------------------------------------------
% Build periodic active row cache (cell-list path)

fprintf('\n[periodic] building solver row cache...\n');
tRow = tic;
rowCache = geom.build_active_row_cache_periodic( ...
    polsys, problem, ewaldParams, scfParams);
timeRow = toc(tRow);

fprintf('[periodic] solver row cache done in %.6f s | nActive = %d\n', ...
    timeRow, numel(problem.activeSites));

%% ------------------------------------------------------------------------
% Build k-space cache

fprintf('\n[periodic] building k-space cache...\n');
tK = tic;
kCache = geom.build_periodic_kspace_cache( ...
    polsys, problem, ewaldParams, scfParams);
timeK = toc(tK);

fprintf('[periodic] k-space cache done in %.6f s\n', timeK);

%% ------------------------------------------------------------------------
% One induced-field apply benchmark on a probe dipole

rng(1);
mu_probe = zeros(polsys.n_sites, 3);
mu_probe(problem.polMask, :) = cfg.probeDipoleScale * randn(nnz(problem.polMask), 3);

dipoleParams = struct();
dipoleParams.target_mask = targetMask;
dipoleParams.source_mask = targetMask;
dipoleParams.use_thole = cfg.use_thole;
dipoleParams.problem = problem;
dipoleParams.realspace_row_cache = rowCache;
dipoleParams.kspace_cache = kCache;
dipoleParams.kspace_mode = cfg.kspace_mode;
dipoleParams.kspace_memory_limit_gb = cfg.kspace_memory_limit_gb;
dipoleParams.k_block_size = cfg.k_block_size;
dipoleParams.verbose = false;

fprintf('\n[periodic] benchmarking one induced-field apply on probe dipoles...\n');

for k = 1:cfg.nWarmup
    thole.induced_field_from_dipoles_thole_periodic( ...
        polsys, mu_probe, ewaldParams, dipoleParams);
end

tApply = zeros(cfg.nRepeat, 1);
for k = 1:cfg.nRepeat
    tic;
    Eprobe = thole.induced_field_from_dipoles_thole_periodic( ...
        polsys, mu_probe, ewaldParams, dipoleParams);
    tApply(k) = toc;
end

fprintf('[periodic] one-apply mean = %.6f s | ||Eprobe||_F = %.16e\n', ...
    mean(tApply), norm(Eprobe, 'fro'));

%% ------------------------------------------------------------------------
% Periodic SOR solve using periodic active row cache

scfParamsSolve = scfParams;
scfParamsSolve.realspace_dipole_row_cache = rowCache;
scfParamsSolve.kspace_cache = kCache;

fprintf('\n[periodic] solving periodic SOR...\n');
tSolve = tic;
[mu, scfOut] = thole.solve_scf_iterative_periodic_sor( ...
    polsys, Eext, ewaldParams, scfParamsSolve);
timeSolve = toc(tSolve);

fprintf('[periodic] SOR done in %.6f s | nIter = %d | relres = %.3e\n', ...
    timeSolve, scfOut.nIter, scfOut.relres);

%% ------------------------------------------------------------------------
% Dipole summary

local_print_dipole_summary('periodic benchmark (cell-list active row)', polsys, mu);

fprintf('\n--- timing summary ---\n');
fprintf('  periodic Eext         = %.6f s\n', timeEext);
fprintf('  problem prep          = %.6f s\n', timeProblem);
fprintf('  solver row cache      = %.6f s\n', timeRow);
fprintf('  k-space cache         = %.6f s\n', timeK);
fprintf('  induced-field apply   = %.6f s\n', mean(tApply));
fprintf('  periodic SOR          = %.6f s\n', timeSolve);

fprintf('\nDone.\n');

%% ========================================================================
% local helpers
%% ========================================================================

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