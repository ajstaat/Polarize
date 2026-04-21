%% periodic_loose_sor_timing_test
% Periodic loose-SOR timing test for manual supercell-size sweeps.
%
% Intended use:
%   - change cfg.supercellSize by hand
%   - run
%   - record timings / nIter / residual / energy
%
% This script:
%   1) builds a charged-pair real-crystal test case
%   2) computes periodic external field
%   3) prepares the periodic SCF problem
%   4) builds shared periodic matrix-free caches once
%   5) runs loose matrix-free periodic SOR
%   6) reports timing and final energy
%
% Notes:
%   - Uses the validated periodic matrix-free SOR path
%   - Uses shared-cache inputs so the solver itself does not rebuild caches
%   - No direct solve or GMRES here; this is just the cheap timing driver

clear; clc; close all;

fprintf('=== periodic loose-SOR timing test ===\n');

%% ------------------------------------------------------------------------
% User controls

cfg = struct();

cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename   = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

% Change this manually for timing sweeps
cfg.supercellSize = [2 5 1];

cfg.bondScale   = 1.20;
cfg.pairCharges = [+1 -1];

cfg.use_thole = true;
cfg.verbose   = true;
cfg.softening = 0.0;

cfg.periodic = struct();
cfg.periodic.alpha    = 0.30;
cfg.periodic.rcut     = 12.0;
cfg.periodic.kcut     = 3.5;
cfg.periodic.boundary = 'tinfoil';

% Optional chunk size for periodic charge-field k-space evaluation
cfg.periodic.k_block_size = 2048;

% Loose SOR settings
cfg.sor = struct();
cfg.sor.tol           = 1e-6;
cfg.sor.maxIter       = 200;
cfg.sor.omega         = 0.45;
cfg.sor.printEvery    = 10;
cfg.sor.residualEvery = 10;
cfg.sor.stopMetric    = 'max_dmu';

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

targetMask = logical(polsys.site_is_polarizable(:));
sourceMask = abs(polsys.site_charge(:)) > 0;

fprintf('\nMask summary:\n');
fprintf('  n target sites (polarizable) = %d\n', nnz(targetMask));
fprintf('  n source sites (charged)     = %d\n', nnz(sourceMask));
fprintf('  total source charge          = %+0.12f\n', sum(polsys.site_charge(sourceMask)));
fprintf('  max |source charge|          = %.12e\n', max(abs(polsys.site_charge(sourceMask))));

%% ------------------------------------------------------------------------
% Shared periodic params

ewaldParams = struct();
ewaldParams.alpha    = cfg.periodic.alpha;
ewaldParams.rcut     = cfg.periodic.rcut;
ewaldParams.kcut     = cfg.periodic.kcut;
ewaldParams.boundary = cfg.periodic.boundary;

fprintf('\n============================================================\n');
fprintf('PERIODIC LOOSE-SOR SETTINGS\n');
fprintf('============================================================\n');
fprintf('  supercell = [%d %d %d]\n', cfg.supercellSize);
fprintf('  alpha     = %.6f\n', ewaldParams.alpha);
fprintf('  rcut      = %.6f bohr\n', ewaldParams.rcut);
fprintf('  kcut      = %.6f bohr^-1\n', ewaldParams.kcut);
fprintf('  boundary  = %s\n', ewaldParams.boundary);
fprintf('  tol       = %.3e\n', cfg.sor.tol);
fprintf('  maxIter   = %d\n', cfg.sor.maxIter);
fprintf('  omega     = %.3f\n', cfg.sor.omega);
fprintf('  stopMetric= %s\n', cfg.sor.stopMetric);

%% ------------------------------------------------------------------------
% Build periodic external field

params_per = struct();
params_per.use_thole = cfg.use_thole;
params_per.field = struct();
params_per.field.mode = 'nonperiodic';
params_per.field.exclude_self = true;
params_per.field.use_thole_damping = cfg.use_thole;
params_per.field.target_mask = targetMask;
params_per.field.source_mask = sourceMask;
params_per.field.ewald = ewaldParams;
params_per.field.k_block_size = cfg.periodic.k_block_size;

fprintf('\n[periodic] building external field...\n');
tField = tic;
Eext = calc.compute_external_field(polsys, params_per);
timeField = toc(tField);
fprintf('[periodic] external field done in %.6f s | ||Eext||_F = %.16e\n', ...
    timeField, norm(Eext, 'fro'));

%% ------------------------------------------------------------------------
% Prepare SCF problem

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
scfParams.kspace_mode = 'auto';
scfParams.kspace_memory_limit_gb = 8;
scfParams.k_block_size = 2048;

fprintf('[periodic] preparing SCF problem...\n');
tProblem = tic;
problem = thole.prepare_scf_problem(polsys, Eext, scfParams);
timeProblem = toc(tProblem);
fprintf('[periodic] problem prep done in %.6f s | nPolSites = %d\n', ...
    timeProblem, problem.nPolSites);

%% ------------------------------------------------------------------------
% Build shared matrix-free caches once, with separate timings

fprintf('[periodic] building shared matrix-free caches...\n');

cacheParams = struct();
cacheParams.use_thole = cfg.use_thole;
cacheParams.verbose = true;

tRealCache = tic;
realCache_shared = geom.build_periodic_realspace_cache( ...
    polsys, problem, ewaldParams, cacheParams);
timeRealCache = toc(tRealCache);
fprintf('  real-space cache done in %.6f s | nReal=%d\n', ...
    timeRealCache, realCache_shared.nInteractions);

tRowCache = tic;
rowCache_shared = geom.build_periodic_realspace_row_cache( ...
    polsys, problem, realCache_shared);
timeRowCache = toc(tRowCache);
fprintf('  row cache done in %.6f s | nDirected=%d\n', ...
    timeRowCache, rowCache_shared.nInteractions);

tKCache = tic;
kCache_shared = geom.build_periodic_kspace_cache( ...
    polsys, problem, ewaldParams, scfParams);
timeKCache = toc(tKCache);
fprintf('  k-space cache done in %.6f s | nK=%d\n', ...
    timeKCache, kCache_shared.num_kvec);

timeCache = timeRealCache + timeRowCache + timeKCache;
fprintf('[periodic] shared cache build total = %.6f s\n', timeCache);

%% ------------------------------------------------------------------------
% Loose matrix-free periodic SOR solve

fprintf('\n============================================================\n');
fprintf('PERIODIC MATRIX-FREE LOOSE SOR SOLVE\n');
fprintf('============================================================\n');

scfParams.realspace_cache = realCache_shared;
scfParams.realspace_row_cache = rowCache_shared;
scfParams.kspace_cache = kCache_shared;

fprintf('[sor] solving periodic matrix-free SOR...\n');
tSOR = tic;
[mu_sor, scf_sor] = thole.solve_scf_iterative_periodic_sor( ...
    polsys, Eext, ewaldParams, scfParams);
timeSOR = toc(tSOR);
fprintf('[sor] solve done in %.6f s | ||mu||_2 = %.16e\n', ...
    timeSOR, norm(mu_sor));
fprintf('  converged = %d | nIter = %d | relres(iterative) = %.16e\n', ...
    scf_sor.converged, scf_sor.nIter, scf_sor.relres);

%% ------------------------------------------------------------------------
% Post-solve energy / residual check without assembling Tpol

fprintf('[sor] computing energy / physical residual...\n');
tPost = tic;

dipoleParams = struct();
dipoleParams.use_thole = cfg.use_thole;
dipoleParams.problem = problem;
dipoleParams.target_mask = problem.polMask;
dipoleParams.source_mask = problem.polMask;
dipoleParams.realspace_cache = realCache_shared;
dipoleParams.kspace_cache = kCache_shared;

Edip = thole.induced_field_from_dipoles_thole_periodic( ...
    polsys, mu_sor, ewaldParams, dipoleParams);

rhs = zeros(problem.nSites, 3);
rhs(problem.polMask, :) = problem.alpha(problem.polMask) .* ...
    (Eext(problem.polMask, :) + Edip(problem.polMask, :));

r = rhs(problem.polMask, :) - mu_sor(problem.polMask, :);
bn = norm(rhs(problem.polMask, :), 'fro');
if bn == 0
    bn = 1.0;
end
relres_phys = norm(r, 'fro') / bn;

polMask = problem.polMask;
mu_pol = mu_sor(polMask, :);
Eext_pol = Eext(polMask, :);
Edip_pol = Edip(polMask, :);

E_total = -0.5 * sum(sum(mu_pol .* Eext_pol));

timePost = toc(tPost);

fprintf('[sor] postprocessing done in %.6f s\n', timePost);
fprintf('  relres(physical) = %.16e\n', relres_phys);
fprintf('  total energy     = %+0.12e Ha\n', E_total);

%% ------------------------------------------------------------------------
% Timing summary

fprintf('\n============================================================\n');
fprintf('TIMING SUMMARY\n');
fprintf('============================================================\n');
fprintf('  field_time_s      = %.6f\n', timeField);
fprintf('  problem_time_s    = %.6f\n', timeProblem);
fprintf('  real_cache_s      = %.6f\n', timeRealCache);
fprintf('  row_cache_s       = %.6f\n', timeRowCache);
fprintf('  k_cache_s         = %.6f\n', timeKCache);
fprintf('  shared_cache_s    = %.6f\n', timeCache);
fprintf('  sor_solve_s       = %.6f\n', timeSOR);
fprintf('  postprocess_s     = %.6f\n', timePost);

q = polsys.site_charge(:);
r = polsys.site_pos;

p_src = sum(q .* r, 1);          % e*bohr
p_ind = sum(mu_sor, 1);          % e*bohr
p_tot = p_src + p_ind;           % e*bohr

convD = 2.541746;

fprintf('Bare source dipole  (e*bohr): [% .6f % .6f % .6f]\n', p_src);
fprintf('Bare source |p|     = %.6f e*bohr = %.6f D\n', norm(p_src), norm(p_src)*convD);

fprintf('Induced dipole sum  (e*bohr): [% .6f % .6f % .6f]\n', p_ind);
fprintf('Induced    |p|      = %.6f e*bohr = %.6f D\n', norm(p_ind), norm(p_ind)*convD);

fprintf('Total cell dipole   (e*bohr): [% .6f % .6f % .6f]\n', p_tot);
fprintf('Total      |p|      = %.6f e*bohr = %.6f D\n', norm(p_tot), norm(p_tot)*convD);

polMask = logical(polsys.site_is_polarizable(:));

mu_pol   = mu_sor(polMask, :);
Eext_pol = Eext(polMask, :);
Edip_pol = Edip(polMask, :);
Eloc_pol = Eext_pol + Edip_pol;

fprintf('sum |Eext_i|^2   = %.16e\n', sum(sum(Eext_pol.^2)));
fprintf('sum |Edip_i|^2   = %.16e\n', sum(sum(Edip_pol.^2)));
fprintf('sum |Eloc_i|^2   = %.16e\n', sum(sum(Eloc_pol.^2)));

fprintf('sum mu_i·Eext_i  = %.16e\n', sum(sum(mu_pol .* Eext_pol)));
fprintf('sum mu_i·Edip_i  = %.16e\n', sum(sum(mu_pol .* Edip_pol)));
fprintf('sum mu_i·Eloc_i  = %.16e\n', sum(sum(mu_pol .* Eloc_pol)));

fprintf('sum |mu_i|^2     = %.16e\n', sum(sum(mu_pol.^2)));
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