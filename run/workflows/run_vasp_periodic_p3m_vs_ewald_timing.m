%% run_vasp_periodic_p3m_vs_ewald_timing
%
% Real VASP [2 5 1] periodic timing workflow:
%
%   VASP/CONTCAR
%   -> crystal template
%   -> periodic supercell system
%   -> centered charged molecular pair
%   -> periodic polarization system
%   -> Ewald external field + periodic_ewald GMRES solve
%   -> P3M external field   + periodic_p3m GMRES solve
%   -> timing / apply / energy comparison
%
% This is deliberately modeled after run_vasp_periodic_polarization_workflow.
%
% Periodic convention:
%
%   Polarize uses direct lattice vectors as ROWS:
%
%       cart = frac * H
%
% Public workflow code should not manually transpose the lattice.

clear; clc; close all;

fprintf('\n============================================================\n');
fprintf('Real VASP periodic P3M vs Ewald timing workflow\n');
fprintf('============================================================\n');

HARTREE_TO_EV = 27.211386245988;

%% ------------------------------------------------------------------------
% User controls
% -------------------------------------------------------------------------

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

% -------------------------------------------------------------------------
% Periodic external fixed-charge field: Ewald reference.
% -------------------------------------------------------------------------

cfg.field = struct();
cfg.field.mode = 'periodic';
cfg.field.exclude_self = true;
cfg.field.use_thole_damping = true;
cfg.field.kspace_mode = 'auto';      % 'auto' | 'full' | 'blocked'
cfg.field.k_block_size = 2048;
cfg.field.kspace_memory_limit_gb = 8;
cfg.field.real_only = false;
cfg.field.verbose = false;

cfg.field.ewald = struct();
cfg.field.ewald.alpha = 0.30;
cfg.field.ewald.rcut = 18.0;
cfg.field.ewald.kcut = 2.5;
cfg.field.ewald.boundary = 'tinfoil';

% -------------------------------------------------------------------------
% P3M controls.
% -------------------------------------------------------------------------

cfg.p3m = struct();

% First-pass mesh. Tune this after seeing timing/error.
cfg.p3m.mesh_size = [32 48 32];
cfg.p3m.assignment_order = 4;

cfg.p3m.derivative_mode = 'spectral';
cfg.p3m.influence_mode = 'ewald';
cfg.p3m.fd_stencil = 'central2';

cfg.p3m.deconvolve_assignment = true;
cfg.p3m.deconvolution_floor = 1e-8;
cfg.p3m.alias_range = 2;

% -------------------------------------------------------------------------
% SCF / solver controls.
% -------------------------------------------------------------------------

cfg.scf = struct();
cfg.scf.tol = 1e-8;
cfg.scf.maxIter = 500;
cfg.scf.mixing = 0.6;
cfg.scf.omega = 1.0;
cfg.scf.verbose = false;

% This timing workflow compares GMRES/global apply, because that is the
% validated P3M path. SOR/P3M fast path can be tested separately later.
cfg.solver = struct();
cfg.solver.method = 'gmres';
cfg.solver.gmres_restart = [];
cfg.solver.compute_residual = true;
cfg.solver.jacobi_mixing = 0.6;
cfg.solver.sor_omega = 0.97;
cfg.solver.stop_metric = 'max_dmu';
cfg.solver.sor_residual_every = 25;

% -------------------------------------------------------------------------
% Operator factory controls.
% -------------------------------------------------------------------------

cfg.operator = struct();

cfg.operator.backend = 'auto';
cfg.operator.use_thole = true;
cfg.operator.softening = 0.0;

cfg.operator.alpha = cfg.field.ewald.alpha;
cfg.operator.rcut = cfg.field.ewald.rcut;
cfg.operator.kcut = cfg.field.ewald.kcut;
cfg.operator.boundary = cfg.field.ewald.boundary;

cfg.operator.kspace_mode = cfg.field.kspace_mode;
cfg.operator.k_block_size = cfg.field.k_block_size;
cfg.operator.kspace_memory_limit_gb = cfg.field.kspace_memory_limit_gb;

cfg.operator.use_mex = true;
cfg.operator.use_mex_kspace = true;
cfg.operator.profile = false;
cfg.operator.verbose = true;

% Operator-apply timing.
cfg.timing = struct();
cfg.timing.apply_repeats = 5;

%% ------------------------------------------------------------------------
% Check file and print controls
% -------------------------------------------------------------------------

if ~isfile(cfg.filename)
    error('VASP file not found:\n  %s\nEdit cfg.filename at the top of this script.', cfg.filename);
end

fprintf('\nInput structure:\n  %s\n', cfg.filename);

fprintf('\nRun controls:\n');
fprintf('  supercell              = [%d %d %d]\n', cfg.supercellSize);
fprintf('  relation/shell         = %s / %d\n', cfg.relation, cfg.shell);
fprintf('  stackAxis/direction    = %s / %s\n', cfg.stackAxis, cfg.direction);
fprintf('  charges                = [%+.3f %+.3f]\n', cfg.pairCharges);
fprintf('  solver method          = %s\n', cfg.solver.method);
fprintf('  Ewald alpha            = %.6g\n', cfg.operator.alpha);
fprintf('  Ewald rcut             = %.6g bohr\n', cfg.operator.rcut);
fprintf('  Ewald kcut             = %.6g bohr^-1\n', cfg.operator.kcut);
fprintf('  Ewald boundary         = %s\n', cfg.operator.boundary);
fprintf('  Ewald kspace mode      = %s\n', cfg.operator.kspace_mode);
fprintf('  P3M mesh               = [%d %d %d]\n', cfg.p3m.mesh_size);
fprintf('  P3M assignment order   = %d\n', cfg.p3m.assignment_order);
fprintf('  field Thole damping    = %d\n', cfg.field.use_thole_damping);
fprintf('  operator Thole damping = %d\n', cfg.operator.use_thole);
fprintf('  operator use_mex       = %d\n', cfg.operator.use_mex);

%% ------------------------------------------------------------------------
% 1. Import crystal template
% -------------------------------------------------------------------------

fprintf('\n[1] Importing crystal template...\n');

tImport = tic;
crystal = io.import_contcar_as_crystal(cfg.filename, ...
    'BondScale', cfg.bondScale, ...
    'SortMolecules', false);
importTime = toc(tImport);

fprintf('  import time = %.6f s\n', importTime);
fprintf('  nSites      = %d\n', crystal.nSites);
fprintf('  nBaseMols   = %d\n', crystal.nBaseMols);
fprintf('  units       = %s\n', crystal.units.length);

%% ------------------------------------------------------------------------
% 2. Define model
% -------------------------------------------------------------------------

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
    'C_deg3', 1.334, ...
    'C_deg4', 1.750, ...
    'N', 1.073, ...
    'O', 0.837);

fprintf('  alpha_units = %s\n', model.alpha_units);
fprintf('  thole_a     = %.3f\n', model.thole_a);

%% ------------------------------------------------------------------------
% 3. Build system
% -------------------------------------------------------------------------

fprintf('\n[3] Building crystal system...\n');

buildOpts = struct();
buildOpts.supercell_size = cfg.supercellSize;
buildOpts.bondScale = cfg.bondScale;
buildOpts.verbose = cfg.verbose;

tBuild = tic;
sys0 = builder.make_crystal_system(crystal, model, buildOpts);
buildTime = toc(tBuild);

fprintf('\nBuilt system:\n');
fprintf('  build time       = %.6f s\n', buildTime);
fprintf('  supercell_size   = [%d %d %d]\n', sys0.supercell_size);
fprintf('  n_sites          = %d\n', sys0.n_sites);
fprintf('  n_molecules      = %d\n', numel(sys0.molecule_table.molecule_id));
fprintf('  n_complete       = %d\n', nnz(sys0.molecule_table.is_complete_in_display));
fprintf('  length unit      = %s\n', sys0.units.length);
fprintf('  alpha unit       = %s\n', sys0.units.alpha);

completeIDs = builder.complete_molecule_ids(sys0);
if isempty(completeIDs)
    error('No complete molecules found in the displayed supercell.');
end

lat0 = geom.get_lattice(sys0);
Lmin0 = geom.shortest_lattice_translation(lat0.H);

fprintf('  lattice volume   = %.8e bohr^3\n', lat0.volume);
fprintf('  ||H*G - 2piI||   = %.3e\n', norm(lat0.H * lat0.G - 2*pi*eye(3), 'fro'));
fprintf('  Lmin             = %.8f bohr\n', Lmin0);
fprintf('  Lmin/2           = %.8f bohr\n', 0.5 * Lmin0);

if ~(cfg.field.ewald.rcut < 0.5 * Lmin0)
    error(['Periodic Ewald real-space cache requires rcut < Lmin/2.\n' ...
           '  rcut   = %.8f bohr\n' ...
           '  Lmin/2 = %.8f bohr\n' ...
           'Increase the supercell or reduce cfg.field.ewald.rcut.'], ...
           cfg.field.ewald.rcut, 0.5 * Lmin0);
end

%% ------------------------------------------------------------------------
% 4. Select centered neighbor pair
% -------------------------------------------------------------------------

fprintf('\n[4] Selecting centered %s shell %d pair...\n', cfg.relation, cfg.shell);

tSelect = tic;
selection = builder.select_centered_neighbor_pair(sys0, ...
    'Relation', cfg.relation, ...
    'Shell', cfg.shell, ...
    'StackAxis', cfg.stackAxis, ...
    'Direction', cfg.direction, ...
    'Verbose', cfg.verbose);
selectTime = toc(tSelect);

pairVector0 = selection.neighbor_com - selection.reference_com;
pairDistance0 = norm(pairVector0);

fprintf('  pair selection time = %.6f s\n', selectTime);
fprintf('\nSelected pair before charging:\n');
fprintf('  relation       = %s\n', selection.relation);
fprintf('  reference ID   = %d\n', selection.reference_mol_id);
fprintf('  neighbor ID    = %d\n', selection.neighbor_mol_id);
fprintf('  pair distance  = %.6f bohr\n', pairDistance0);
fprintf('  midpoint dist  = %.6f bohr\n', selection.midpoint_distance);
fprintf('  pair midpoint / bohr:\n');
disp(selection.pair_midpoint);

%% ------------------------------------------------------------------------
% 5. Apply uniform charges
% -------------------------------------------------------------------------

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
fprintf('  selected source charge       = %+ .8f e\n', ...
    sum(sys.site_charge(idxRef)) + sum(sys.site_charge(idxNbr)));
fprintf('  total system charge          = %+ .8f e\n', sum(sys.site_charge));

chargedMask = ismember(sys.site_mol_id, molIDs);

fprintf('\nActive/polarizable summary:\n');
fprintf('  active sites                 = %d\n', nnz(sys.site_is_active));
fprintf('  charged sites                = %d\n', nnz(chargedMask));
fprintf('  charged sites polarizable?   = %d\n', any(sys.site_is_polarizable(chargedMask)));
fprintf('  remaining polarizable sites  = %d\n', nnz(sys.site_is_polarizable));

fprintf('\nCharged-source alpha diagnostic:\n');
fprintf('  nnz charged site_alpha       = %d / %d\n', ...
    nnz(sys.site_alpha(chargedMask) > 0), nnz(chargedMask));
fprintf('  min charged site_alpha       = %.12e\n', min(sys.site_alpha(chargedMask)));
fprintf('  max charged site_alpha       = %.12e\n', max(sys.site_alpha(chargedMask)));

if any(sys.site_is_polarizable(chargedMask))
    error('Charged molecule sites should not remain polarizable.');
end

if abs(sum(sys.site_charge(chargedMask))) > 1e-10
    error('Selected charged pair must be neutral for periodic fixed-charge field.');
end

io.assert_atomic_units(sys);

%% ------------------------------------------------------------------------
% 6. Extract periodic polarization system
% -------------------------------------------------------------------------

fprintf('\n[6] Extracting periodic polarization system...\n');

polsys = builder.extract_polarization_system(sys, struct('mode', 'periodic'));
io.assert_atomic_units(polsys);

fprintf('  polsys.n_sites     = %d\n', polsys.n_sites);
fprintf('  polsys.is_periodic = %d\n', polsys.is_periodic);
fprintf('  polarizable sites  = %d\n', nnz(polsys.site_is_polarizable));
fprintf('  charged sites      = %d\n', nnz(abs(polsys.site_charge) > 0));

lat = geom.get_lattice(polsys);
Lmin = geom.shortest_lattice_translation(lat.H);

fprintf('  lattice volume     = %.8e bohr^3\n', lat.volume);
fprintf('  ||H*G - 2piI||     = %.3e\n', norm(lat.H * lat.G - 2*pi*eye(3), 'fro'));
fprintf('  Lmin               = %.8f bohr\n', Lmin);
fprintf('  Lmin/2             = %.8f bohr\n', 0.5 * Lmin);

if ~(cfg.field.ewald.rcut < 0.5 * Lmin)
    error(['Periodic real-space cache requires rcut < Lmin/2 after extraction.\n' ...
           '  rcut   = %.8f bohr\n' ...
           '  Lmin/2 = %.8f bohr'], ...
           cfg.field.ewald.rcut, 0.5 * Lmin);
end

targetMask = logical(polsys.site_is_polarizable(:));
sourceMask = abs(polsys.site_charge(:)) > 0;

if abs(sum(polsys.site_charge(sourceMask))) > 1e-10
    error('Periodic selected sources are not neutral.');
end

%% ------------------------------------------------------------------------
% 7. Ewald external field
% -------------------------------------------------------------------------

fprintf('\n[7] Computing periodic Ewald external field...\n');

fieldParamsEwald = struct();
fieldParamsEwald.use_thole = cfg.field.use_thole_damping;
fieldParamsEwald.field = struct();
fieldParamsEwald.field.mode = cfg.field.mode;
fieldParamsEwald.field.exclude_self = cfg.field.exclude_self;
fieldParamsEwald.field.use_thole_damping = cfg.field.use_thole_damping;
fieldParamsEwald.field.target_mask = targetMask;
fieldParamsEwald.field.source_mask = sourceMask;
fieldParamsEwald.field.real_only = cfg.field.real_only;
fieldParamsEwald.field.kspace_mode = cfg.field.kspace_mode;
fieldParamsEwald.field.k_block_size = cfg.field.k_block_size;
fieldParamsEwald.field.kspace_memory_limit_gb = cfg.field.kspace_memory_limit_gb;
fieldParamsEwald.field.verbose = cfg.field.verbose;
fieldParamsEwald.field.ewald = cfg.field.ewald;

tFieldEwald = tic;
EextEwald = calc.compute_external_field(polsys, fieldParamsEwald);
fieldTimeEwald = toc(tFieldEwald);

fieldDirect = fieldParamsEwald.field;
fieldDirect = rmfield(fieldDirect, 'mode');
[~, fieldPartsEwald] = thole.induced_field_from_charges_periodic(polsys, fieldDirect);

fprintf('  Ewald Eext computed in %.6f s\n', fieldTimeEwald);
fprintf('  ||Eext||_F                 = %.12e\n', norm(EextEwald, 'fro'));
fprintf('  ||Eext polarizable||_F     = %.12e\n', norm(EextEwald(targetMask, :), 'fro'));
fprintf('  ||Ereal||_F                = %.12e\n', norm(fieldPartsEwald.real, 'fro'));
fprintf('  ||Erecip||_F               = %.12e\n', norm(fieldPartsEwald.recip, 'fro'));
fprintf('  ||Esurf||_F                = %.12e\n', norm(fieldPartsEwald.surf, 'fro'));
fprintf('  field nK                   = %d\n', fieldPartsEwald.nK);
fprintf('  field real entries         = %d\n', fieldPartsEwald.nRealEntries);
fprintf('  field storage mode         = %s\n', fieldPartsEwald.storage_mode);
fprintf('  selected source charge     = %+ .12e\n', fieldPartsEwald.qtot);

if isfield(fieldPartsEwald, 'realCache') && ...
        isfield(fieldPartsEwald.realCache, 'B_ewald') && ...
        isfield(fieldPartsEwald.realCache, 'thole_delta')
    fprintf('  ||field B_ewald||_2        = %.12e\n', norm(fieldPartsEwald.realCache.B_ewald));
    fprintf('  ||field thole_delta||_2    = %.12e\n', norm(fieldPartsEwald.realCache.thole_delta));
    fprintf('  thole/ewald coeff ratio    = %.12e\n', ...
        norm(fieldPartsEwald.realCache.thole_delta) / max(norm(fieldPartsEwald.realCache.B_ewald), eps));
end

%% ------------------------------------------------------------------------
% 8. P3M external field
% -------------------------------------------------------------------------

fprintf('\n[8] Computing periodic P3M external field...\n');

fieldParamsP3M = struct();
fieldParamsP3M.ewald = cfg.field.ewald;
fieldParamsP3M.mesh_size = cfg.p3m.mesh_size;
fieldParamsP3M.assignment_order = cfg.p3m.assignment_order;
fieldParamsP3M.target_mask = targetMask;
fieldParamsP3M.source_mask = sourceMask;
fieldParamsP3M.exclude_self = cfg.field.exclude_self;
fieldParamsP3M.use_thole_damping = cfg.field.use_thole_damping;
fieldParamsP3M.realspace_backend = 'thole_periodic_real';
fieldParamsP3M.derivative_mode = cfg.p3m.derivative_mode;
fieldParamsP3M.influence_mode = cfg.p3m.influence_mode;
fieldParamsP3M.fd_stencil = cfg.p3m.fd_stencil;
fieldParamsP3M.deconvolve_assignment = cfg.p3m.deconvolve_assignment;
fieldParamsP3M.deconvolution_floor = cfg.p3m.deconvolution_floor;
fieldParamsP3M.alias_range = cfg.p3m.alias_range;
fieldParamsP3M.verbose = cfg.field.verbose;

tFieldP3M = tic;
[EextP3M, fieldPartsP3M] = p3m.compute_external_field_charges(polsys, fieldParamsP3M);
fieldTimeP3M = toc(tFieldP3M);

EextDiff = EextP3M(targetMask, :) - EextEwald(targetMask, :);
EextRelDiff = norm(EextDiff, 'fro') / max(norm(EextEwald(targetMask, :), 'fro'), eps);
EextCos = local_cosine(EextP3M(targetMask, :), EextEwald(targetMask, :));

fprintf('  P3M Eext computed in %.6f s\n', fieldTimeP3M);
fprintf('  ||Eext||_F                 = %.12e\n', norm(EextP3M, 'fro'));
fprintf('  ||Eext polarizable||_F     = %.12e\n', norm(EextP3M(targetMask, :), 'fro'));
fprintf('  ||Ereal||_F                = %.12e\n', norm(fieldPartsP3M.real, 'fro'));
fprintf('  ||Erecip||_F               = %.12e\n', norm(fieldPartsP3M.recip, 'fro'));
fprintf('  ||Esurf||_F                = %.12e\n', norm(fieldPartsP3M.surf, 'fro'));
fprintf('  field nK mesh              = %d\n', fieldPartsP3M.nK);
fprintf('  assigned rho total         = %+ .12e\n', fieldPartsP3M.rho_total);
fprintf('  Eext P3M-Ewald rel diff    = %.12e\n', EextRelDiff);
fprintf('  Eext P3M/Ewald cosine      = %.12f\n', EextCos);

%% ------------------------------------------------------------------------
% 9. Prepare SCF problems
% -------------------------------------------------------------------------

fprintf('\n[9] Preparing SCF problems...\n');

problemEwald = thole.prepare_scf_problem(polsys, EextEwald, cfg.scf);
problemP3M = thole.prepare_scf_problem(polsys, EextP3M, cfg.scf);

fprintf('  nPolSites                  = %d\n', problemEwald.nPolSites);
fprintf('  active vector length       = %d\n', numel(problemEwald.Eext_pol_vec));
fprintf('  ||Eext_ewald_pol_vec||     = %.12e\n', norm(problemEwald.Eext_pol_vec));
fprintf('  ||Eext_p3m_pol_vec||       = %.12e\n', norm(problemP3M.Eext_pol_vec));
fprintf('  rel Eext_pol_vec diff      = %.12e\n', ...
    norm(problemP3M.Eext_pol_vec - problemEwald.Eext_pol_vec) / max(norm(problemEwald.Eext_pol_vec), eps));

%% ------------------------------------------------------------------------
% 10. Build periodic Ewald operator
% -------------------------------------------------------------------------

fprintf('\n[10] Building periodic Ewald solver operator...\n');

tOpEwald = tic;
opEwald = thole.make_polarization_operator(polsys, problemEwald, ...
    'Mode', 'periodic_ewald', ...
    'Solver', cfg.solver.method, ...
    'Backend', cfg.operator.backend, ...
    'UseThole', cfg.operator.use_thole, ...
    'Softening', cfg.operator.softening, ...
    'Rcut', cfg.operator.rcut, ...
    'Alpha', cfg.operator.alpha, ...
    'Kcut', cfg.operator.kcut, ...
    'Boundary', cfg.operator.boundary, ...
    'KspaceMode', cfg.operator.kspace_mode, ...
    'KBlockSize', cfg.operator.k_block_size, ...
    'KspaceMemoryLimitGB', cfg.operator.kspace_memory_limit_gb, ...
    'UseMex', cfg.operator.use_mex, ...
    'UseMexKspace', cfg.operator.use_mex_kspace, ...
    'Profile', cfg.operator.profile, ...
    'Verbose', cfg.operator.verbose);
opTimeEwald = toc(tOpEwald);

local_print_operator_summary(opEwald, opTimeEwald, 'Ewald');

%% ------------------------------------------------------------------------
% 11. Build periodic P3M operator
% -------------------------------------------------------------------------

fprintf('\n[11] Building periodic P3M solver operator...\n');

tOpP3M = tic;
opP3M = thole.make_polarization_operator(polsys, problemP3M, ...
    'Mode', 'periodic_p3m', ...
    'Solver', cfg.solver.method, ...
    'Backend', cfg.operator.backend, ...
    'UseThole', cfg.operator.use_thole, ...
    'Softening', cfg.operator.softening, ...
    'Rcut', cfg.operator.rcut, ...
    'Alpha', cfg.operator.alpha, ...
    'Kcut', cfg.operator.kcut, ...
    'Boundary', cfg.operator.boundary, ...
    'MeshSize', cfg.p3m.mesh_size, ...
    'AssignmentOrder', cfg.p3m.assignment_order, ...
    'DerivativeMode', cfg.p3m.derivative_mode, ...
    'InfluenceMode', cfg.p3m.influence_mode, ...
    'FDStencil', cfg.p3m.fd_stencil, ...
    'DeconvolveAssignment', cfg.p3m.deconvolve_assignment, ...
    'DeconvolutionFloor', cfg.p3m.deconvolution_floor, ...
    'AliasRange', cfg.p3m.alias_range, ...
    'UseMex', cfg.operator.use_mex, ...
    'Profile', cfg.operator.profile, ...
    'Verbose', cfg.operator.verbose);
opTimeP3M = toc(tOpP3M);

local_print_operator_summary(opP3M, opTimeP3M, 'P3M');

if isfield(opP3M, 'p3m_cache') && isfield(opP3M.p3m_cache, 'estimated_mesh_gb')
    fprintf('  P3M estimated mesh storage = %.6f GB\n', opP3M.p3m_cache.estimated_mesh_gb);
end

%% ------------------------------------------------------------------------
% 12. Operator apply timing
% -------------------------------------------------------------------------

fprintf('\n[12] Timing operator apply...\n');

rng(123);
xTest = randn(numel(problemEwald.Eext_pol_vec), 1);

applyEwald = local_time_apply(opEwald, xTest, cfg.timing.apply_repeats);
applyP3M = local_time_apply(opP3M, xTest, cfg.timing.apply_repeats);

TxEwald = opEwald.apply(xTest);
TxP3M = opP3M.apply(xTest);

applyFullRelDiff = norm(TxP3M - TxEwald) / max(norm(TxEwald), eps);
applyCos = local_cosine(TxP3M, TxEwald);

fprintf('  Ewald apply median         = %.6f s\n', applyEwald.median);
fprintf('  P3M apply median           = %.6f s\n', applyP3M.median);
fprintf('  apply speedup Ewald/P3M    = %.6f x\n', applyEwald.median / max(applyP3M.median, eps));
fprintf('  apply full rel diff        = %.12e\n', applyFullRelDiff);
fprintf('  apply cosine               = %.12f\n', applyCos);

%% ------------------------------------------------------------------------
% 13. Solve Ewald SCF
% -------------------------------------------------------------------------

fprintf('\n[13] Solving Ewald SCF with %s...\n', cfg.solver.method);

tSolveEwald = tic;
[muEwald, solveInfoEwald] = local_solve_selected(problemEwald, opEwald, cfg);
solveWallEwald = toc(tSolveEwald);

local_print_solver_summary(solveInfoEwald, solveWallEwald, muEwald);

%% ------------------------------------------------------------------------
% 14. Solve P3M SCF
% -------------------------------------------------------------------------

fprintf('\n[14] Solving P3M SCF with %s...\n', cfg.solver.method);

tSolveP3M = tic;
[muP3M, solveInfoP3M] = local_solve_selected(problemP3M, opP3M, cfg);
solveWallP3M = toc(tSolveP3M);

local_print_solver_summary(solveInfoP3M, solveWallP3M, muP3M);

%% ------------------------------------------------------------------------
% 15. Energies
% -------------------------------------------------------------------------

fprintf('\n[15] Computing active-space polarization energies...\n');

energyEwald = calc.compute_total_energy_active_space( ...
    polsys, problemEwald, muEwald, EextEwald, opEwald);

energyP3M = calc.compute_total_energy_active_space( ...
    polsys, problemP3M, muP3M, EextP3M, opP3M);

fprintf('\nEwald energy breakdown:\n');
local_print_energy(energyEwald, HARTREE_TO_EV);

fprintf('\nP3M energy breakdown:\n');
local_print_energy(energyP3M, HARTREE_TO_EV);

muRelDiff = norm(muP3M - muEwald, 'fro') / max(norm(muEwald, 'fro'), eps);
muCos = local_cosine(muP3M, muEwald);

fprintf('\nEwald/P3M solution comparison:\n');
fprintf('  ||mu Ewald||_F             = %.12e\n', norm(muEwald, 'fro'));
fprintf('  ||mu P3M||_F               = %.12e\n', norm(muP3M, 'fro'));
fprintf('  rel mu diff                = %.12e\n', muRelDiff);
fprintf('  mu cosine                  = %.12f\n', muCos);
fprintf('  energy diff total          = %+ .12e Ha (%+ .8f eV)\n', ...
    energyP3M.total - energyEwald.total, ...
    (energyP3M.total - energyEwald.total) * HARTREE_TO_EV);
fprintf('  energy diff stationary     = %+ .12e Ha (%+ .8f eV)\n', ...
    energyP3M.total_stationary - energyEwald.total_stationary, ...
    (energyP3M.total_stationary - energyEwald.total_stationary) * HARTREE_TO_EV);

%% ------------------------------------------------------------------------
% 16. Timing summary table
% -------------------------------------------------------------------------

totalEwald = fieldTimeEwald + opTimeEwald + solveWallEwald;
totalP3M = fieldTimeP3M + opTimeP3M + solveWallP3M;

summary = table();
summary.method = ["ewald"; "p3m"];
summary.field_time_s = [fieldTimeEwald; fieldTimeP3M];
summary.operator_build_time_s = [opTimeEwald; opTimeP3M];
summary.apply_median_s = [applyEwald.median; applyP3M.median];
summary.solve_time_s = [solveWallEwald; solveWallP3M];
summary.total_field_op_solve_s = [totalEwald; totalP3M];
summary.energy_total_eV = [energyEwald.total; energyP3M.total] * HARTREE_TO_EV;
summary.energy_stationary_eV = [energyEwald.total_stationary; energyP3M.total_stationary] * HARTREE_TO_EV;
summary.solver_relres = [local_get_numeric_field(solveInfoEwald, 'relres', NaN); ...
                         local_get_numeric_field(solveInfoP3M, 'relres', NaN)];
summary.converged = [local_get_logical_field(solveInfoEwald, 'converged', false); ...
                     local_get_logical_field(solveInfoP3M, 'converged', false)];

fprintf('\n============================================================\n');
fprintf('Periodic Ewald vs P3M timing summary\n');
fprintf('============================================================\n');
disp(summary);

fprintf('\nSpeedups / differences:\n');
fprintf('  field speedup Ewald/P3M       = %.6f x\n', fieldTimeEwald / max(fieldTimeP3M, eps));
fprintf('  op build speedup Ewald/P3M    = %.6f x\n', opTimeEwald / max(opTimeP3M, eps));
fprintf('  apply speedup Ewald/P3M       = %.6f x\n', applyEwald.median / max(applyP3M.median, eps));
fprintf('  solve speedup Ewald/P3M       = %.6f x\n', solveWallEwald / max(solveWallP3M, eps));
fprintf('  total speedup Ewald/P3M       = %.6f x\n', totalEwald / max(totalP3M, eps));
fprintf('  Eext P3M/Ewald rel diff       = %.12e\n', EextRelDiff);
fprintf('  apply P3M/Ewald full rel diff = %.12e\n', applyFullRelDiff);
fprintf('  mu P3M/Ewald rel diff         = %.12e\n', muRelDiff);
fprintf('  Epol P3M-Ewald                = %+ .8f eV\n', ...
    (energyP3M.total - energyEwald.total) * HARTREE_TO_EV);

fprintf('\nRun summary:\n');
fprintf('  file             = %s\n', cfg.filename);
fprintf('  supercell        = [%d %d %d]\n', cfg.supercellSize);
fprintf('  relation/shell   = %s / %d\n', cfg.relation, cfg.shell);
fprintf('  ref/nbr          = %d / %d\n', selection.reference_mol_id, selection.neighbor_mol_id);
fprintf('  pair distance    = %.8f bohr\n', pairDistance0);
fprintf('  charges          = [%+.3f %+.3f]\n', cfg.pairCharges);
fprintf('  Ewald backend    = %s\n', opEwald.backend);
fprintf('  P3M backend      = %s\n', opP3M.backend);
fprintf('  P3M mesh         = [%d %d %d]\n', cfg.p3m.mesh_size);
fprintf('  alpha            = %.6g\n', cfg.operator.alpha);
fprintf('  rcut             = %.6g bohr\n', cfg.operator.rcut);
fprintf('  kcut             = %.6g bohr^-1\n', cfg.operator.kcut);
fprintf('  boundary         = %s\n', cfg.operator.boundary);

fprintf('\nWorkflow completed successfully.\n');

%% =========================================================================
% Local helpers
% =========================================================================

function [mu, solveInfo] = local_solve_selected(problem, op, cfg)
switch lower(cfg.solver.method)
    case 'direct'
        solveOpts = struct();
        solveOpts.compute_residual = cfg.solver.compute_residual;
        [mu, solveInfo] = thole.solve_scf_direct(problem, op, solveOpts);

    case 'jacobi'
        solveOpts = struct();
        solveOpts.tol = cfg.scf.tol;
        solveOpts.max_iter = cfg.scf.maxIter;
        solveOpts.mixing = cfg.solver.jacobi_mixing;
        solveOpts.stop_metric = cfg.solver.stop_metric;
        solveOpts.verbose = cfg.scf.verbose;
        [mu, solveInfo] = thole.solve_scf_jacobi(problem, op, solveOpts);

    case 'gmres'
        solveOpts = struct();
        solveOpts.tol = cfg.scf.tol;
        solveOpts.max_iter = cfg.scf.maxIter;
        solveOpts.restart = cfg.solver.gmres_restart;
        solveOpts.verbose = cfg.scf.verbose;
        [mu, solveInfo] = thole.solve_scf_gmres(problem, op, solveOpts);

    case 'sor'
        solveOpts = struct();
        solveOpts.tol = cfg.scf.tol;
        solveOpts.max_iter = cfg.scf.maxIter;
        solveOpts.omega = cfg.solver.sor_omega;
        solveOpts.stop_metric = cfg.solver.stop_metric;
        solveOpts.residual_every = cfg.solver.sor_residual_every;
        solveOpts.verbose = cfg.scf.verbose;
        [mu, solveInfo] = thole.solve_scf_sor(problem, op, solveOpts);

    otherwise
        error('Unsupported cfg.solver.method "%s".', cfg.solver.method);
end
end

function stats = local_time_apply(op, x, nRepeat)
times = zeros(nRepeat, 1);

% Warmup.
op.apply(x);

for kk = 1:nRepeat
    t = tic;
    op.apply(x);
    times(kk) = toc(t);
end

stats = struct();
stats.times = times;
stats.min = min(times);
stats.mean = mean(times);
stats.median = median(times);
stats.max = max(times);
end

function local_print_operator_summary(op, elapsed, label)
fprintf('  %s operator build wrapper time = %.6f s\n', label, elapsed);
fprintf('  %s operator mode               = %s\n', label, op.mode);
fprintf('  %s operator kind               = %s\n', label, op.kind);
fprintf('  %s operator backend            = %s\n', label, op.backend);
fprintf('  %s operator size               = %d x %d\n', label, op.size(1), op.size(2));

if isfield(op, 'capabilities')
    fprintf('  %s capability apply            = %d\n', label, op.capabilities.apply);
    fprintf('  %s capability dense_matrix     = %d\n', label, op.capabilities.dense_matrix);
    fprintf('  %s capability row_update       = %d\n', label, op.capabilities.row_update);
end

if isfield(op, 'info')
    info = op.info;

    if isfield(info, 'assembly_time')
        fprintf('  %s opinfo assembly time        = %.6f s\n', label, info.assembly_time);
    end
    if isfield(info, 'cache_time')
        fprintf('  %s opinfo cache time           = %.6f s\n', label, info.cache_time);
    end
    if isfield(info, 'real_cache_time')
        fprintf('  %s real cache time             = %.6f s\n', label, info.real_cache_time);
    end
    if isfield(info, 'k_cache_time')
        fprintf('  %s k cache time                = %.6f s\n', label, info.k_cache_time);
    end
    if isfield(info, 'p3m_cache_time')
        fprintf('  %s p3m cache time              = %.6f s\n', label, info.p3m_cache_time);
    end
    if isfield(info, 'fill_time')
        fprintf('  %s opinfo fill time            = %.6f s\n', label, info.fill_time);
    end
    if isfield(info, 'nPairBlocks')
        fprintf('  %s pair blocks total           = %d\n', label, info.nPairBlocks);
    end
    if isfield(info, 'nPairBlocksKept')
        fprintf('  %s pair blocks kept            = %d\n', label, info.nPairBlocksKept);
    end
    if isfield(info, 'nPairBlocksSkippedCutoff')
        fprintf('  %s skipped by cutoff           = %d\n', label, info.nPairBlocksSkippedCutoff);
    end
    if isfield(info, 'nEntriesDirected')
        fprintf('  %s directed row entries        = %d\n', label, info.nEntriesDirected);
    end
    if isfield(info, 'nRealEntriesDirected')
        fprintf('  %s real directed entries       = %d\n', label, info.nRealEntriesDirected);
    end
    if isfield(info, 'nK')
        fprintf('  %s reciprocal k-vectors        = %d\n', label, info.nK);
    end
    if isfield(info, 'mesh_size')
        fprintf('  %s mesh size                   = [%d %d %d]\n', ...
            label, info.mesh_size(1), info.mesh_size(2), info.mesh_size(3));
    end
    if isfield(info, 'assignment_order')
        fprintf('  %s assignment order            = %d\n', label, info.assignment_order);
    end
    if isfield(info, 'alpha')
        fprintf('  %s alpha                       = %.6g\n', label, info.alpha);
    end
    if isfield(info, 'rcut')
        fprintf('  %s rcut                        = %.6g bohr\n', label, info.rcut);
    end
    if isfield(info, 'kcut')
        fprintf('  %s kcut                        = %.6g bohr^-1\n', label, info.kcut);
    end
    if isfield(info, 'boundary')
        fprintf('  %s boundary                    = %s\n', label, info.boundary);
    end
end
end

function local_print_solver_summary(info, wallTime, mu)
fprintf('  solve wall time       = %.6f s\n', wallTime);

if isfield(info, 'solve_time')
    fprintf('  solver internal time  = %.6f s\n', info.solve_time);
end
if isfield(info, 'setup_time')
    fprintf('  setup time            = %.6f s\n', info.setup_time);
end
if isfield(info, 'residual_time')
    fprintf('  residual time         = %.6f s\n', info.residual_time);
end
if isfield(info, 'final_residual_time')
    fprintf('  final residual time   = %.6f s\n', info.final_residual_time);
end
if isfield(info, 'method')
    fprintf('  method                = %s\n', info.method);
end
if isfield(info, 'tol') && ~isempty(info.tol)
    fprintf('  tolerance             = %.12e\n', info.tol);
end
if isfield(info, 'stop_metric')
    fprintf('  stop metric           = %s\n', info.stop_metric);
end
if isfield(info, 'iterations')
    fprintf('  iterations            = %d\n', info.iterations);
end
if isfield(info, 'max_iter')
    fprintf('  max iterations        = %d\n', info.max_iter);
end
if isfield(info, 'restart')
    if isempty(info.restart)
        fprintf('  gmres restart         = [] unrestarted\n');
    else
        fprintf('  gmres restart         = %d\n', info.restart);
    end
end
if isfield(info, 'iter')
    fprintf('  gmres iter            = [%s]\n', num2str(info.iter));
end
if isfield(info, 'flag')
    fprintf('  gmres flag            = %d\n', info.flag);
end
if isfield(info, 'gmres_relres') && ~isempty(info.gmres_relres) && isfinite(info.gmres_relres)
    fprintf('  gmres relres          = %.12e\n', info.gmres_relres);
end
if isfield(info, 'converged') && ~isempty(info.converged)
    fprintf('  converged             = %d\n', logical(info.converged));
end
if isfield(info, 'relres') && ~isempty(info.relres) && isfinite(info.relres)
    fprintf('  SCF relres            = %.12e\n', info.relres);
else
    fprintf('  SCF relres            = skipped/NaN\n');
end
if isfield(info, 'max_dmu') && ~isempty(info.max_dmu) && isfinite(info.max_dmu)
    fprintf('  max dmu               = %.12e\n', info.max_dmu);
end
if isfield(info, 'operator_backend')
    fprintf('  solver op backend     = %s\n', info.operator_backend);
end
if isfield(info, 'rowcache_fast_path_type')
    fprintf('  rowcache fast type    = %s\n', info.rowcache_fast_path_type);
end

fprintf('  ||mu||_F              = %.12e\n', norm(mu, 'fro'));
fprintf('  max |mu_i|            = %.12e\n', max(vecnorm(mu, 2, 2)));
end

function local_print_energy(energy, HARTREE_TO_EV)
fprintf('  polarization_self       = %+ .12e Ha (%+ .8f eV)\n', ...
    energy.polarization_self, energy.polarization_self * HARTREE_TO_EV);
fprintf('  external_charge_dipole  = %+ .12e Ha (%+ .8f eV)\n', ...
    energy.external_charge_dipole, energy.external_charge_dipole * HARTREE_TO_EV);
fprintf('  dipole_dipole           = %+ .12e Ha (%+ .8f eV)\n', ...
    energy.dipole_dipole, energy.dipole_dipole * HARTREE_TO_EV);
fprintf('  total                   = %+ .12e Ha (%+ .8f eV)\n', ...
    energy.total, energy.total * HARTREE_TO_EV);
fprintf('  total_stationary        = %+ .12e Ha (%+ .8f eV)\n', ...
    energy.total_stationary, energy.total_stationary * HARTREE_TO_EV);
fprintf('  stationary_consistency  = %+ .12e Ha\n', energy.stationary_consistency);
fprintf('  energy relres           = %.12e\n', energy.relres);
end

function c = local_cosine(A, B)
a = A(:);
b = B(:);
den = max(norm(a) * norm(b), eps);
c = dot(a, b) / den;
end

function value = local_get_numeric_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end

function value = local_get_logical_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = logical(s.(name));
else
    value = defaultValue;
end
end