%% run_vasp_periodic_polarization_workflow
%
% Real VASP periodic polarization workflow:
%
%   VASP/CONTCAR
%     -> crystal template
%     -> periodic supercell system
%     -> centered charged molecular pair
%     -> periodic polarization system
%     -> periodic Ewald external field from fixed charges
%     -> SCF problem
%     -> solver-aware periodic Ewald polarization operator
%     -> direct / Jacobi / GMRES / SOR solve
%     -> active-space polarization energy
%     -> selected-pair + induced-dipole visualization
%
% Solver/operator contract:
%
%   cfg.solver.method = 'direct'
%       Backend='auto' -> dense_matrix / periodic_paircache_dense
%
%   cfg.solver.method = 'jacobi'
%       Backend='auto' -> matrix_free / periodic_paircache_apply
%
%   cfg.solver.method = 'gmres'
%       Backend='auto' -> matrix_free / periodic_paircache_apply
%
%   cfg.solver.method = 'sor'
%       Backend='auto' -> matrix_free / periodic_rowcache_apply
%
% Periodic convention:
%
%   Polarize uses direct lattice vectors as ROWS:
%
%       cart = frac * H
%
%   The periodic Ewald routines use geom.get_lattice and
%   ewald.enumerate_kvecs_from_lattice. Public workflow code should not
%   manually transpose the lattice.

clear;
clc;
close all;

fprintf('\n============================================================\n');
fprintf('Real VASP periodic polarization workflow\n');
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

% -------------------------------------------------------------------------
% Periodic external fixed-charge field.
%
% The fixed-charge external field uses:
%
%   periodic charge Ewald + local charge-field Thole correction
%
% The induced-dipole operator below uses:
%
%   periodic dipole Ewald + local dipole-field Thole correction
%
% Keeping both enabled is the consistent Thole model path.
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
cfg.field.ewald.kcut = 1.25;
cfg.field.ewald.boundary = 'tinfoil'; % 'tinfoil' | 'vacuum'

% SCF problem defaults.
cfg.scf = struct();
cfg.scf.tol = 1e-8;
cfg.scf.maxIter = 500;
cfg.scf.mixing = 0.6;
cfg.scf.omega = 1.0;
cfg.scf.verbose = false;

% Solver selection:
%   'direct' | 'jacobi' | 'gmres' | 'sor'
cfg.solver = struct();
cfg.solver.method = 'sor';

% Solver-specific overrides.
cfg.solver.compute_residual = true;  % direct only
cfg.solver.jacobi_mixing = 0.6;
cfg.solver.gmres_restart = [];       % [] = unrestarted
cfg.solver.sor_omega = 0.97;

% Iterative fixed-point stopping metric.
%
% Applies to Jacobi and SOR:
%   'relres'  : stop on ||mu - A(E + Tmu)|| / ||A E||
%   'max_dmu' : stop on max per-site ||mu_i(new)-mu_i(old)||
%
% GMRES remains residual-driven by MATLAB gmres.
cfg.solver.stop_metric = 'max_dmu';

% SOR residual diagnostics cadence.
%
% Applies only to SOR when stop_metric='max_dmu'. If stop_metric='relres',
% SOR must compute the residual every iteration because relres is the
% stopping criterion.
%
%   1      : residual every iteration
%   25     : residual on iteration 1, 25, 50, ...
%   0/Inf  : final residual only
cfg.solver.sor_residual_every = 25;

% Operator factory controls.
%
% Public backend values:
%   'auto' | 'dense' | 'matrix_free'
cfg.operator = struct();
cfg.operator.mode = 'periodic_ewald';
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

% Visualization.
cfg.plot = struct();
cfg.plot.showCOM = true;
cfg.plot.showLabels = true;
cfg.plot.showPairLine = true;
cfg.plot.drawBox = true;

cfg.plot.showDipoles = true;
cfg.plot.onlyPolarizableDipoles = true;
cfg.plot.maxArrows = 300;
cfg.plot.dipoleThreshold = 0.0;
cfg.plot.arrowScale = 120.0;      % visual scale only
cfg.plot.arrowLineWidth = 1.0;

HARTREE_TO_EV = 27.211386245988;

%% Check file

if ~isfile(cfg.filename)
    error('VASP file not found:\n  %s\nEdit cfg.filename at the top of this script.', cfg.filename);
end

fprintf('\nInput structure:\n  %s\n', cfg.filename);

fprintf('\nRun controls:\n');
fprintf('  solver method          = %s\n', cfg.solver.method);
fprintf('  operator mode          = %s\n', cfg.operator.mode);
fprintf('  operator backend       = %s\n', cfg.operator.backend);
fprintf('  Ewald alpha            = %.6g\n', cfg.operator.alpha);
fprintf('  Ewald rcut             = %.6g bohr\n', cfg.operator.rcut);
fprintf('  Ewald kcut             = %.6g bohr^-1\n', cfg.operator.kcut);
fprintf('  Ewald boundary         = %s\n', cfg.operator.boundary);
fprintf('  kspace mode            = %s\n', cfg.operator.kspace_mode);
fprintf('  operator use_mex       = %d\n', cfg.operator.use_mex);
fprintf('  field Thole damping    = %d\n', cfg.field.use_thole_damping);
fprintf('  operator Thole damping = %d\n', cfg.operator.use_thole);

%% 1. Import crystal template

fprintf('\n[1] Importing crystal template...\n');

crystal = io.import_contcar_as_crystal(cfg.filename, ...
    'BondScale', cfg.bondScale, ...
    'SortMolecules', false);

fprintf('  nSites    = %d\n', crystal.nSites);
fprintf('  nBaseMols = %d\n', crystal.nBaseMols);
fprintf('  units     = %s\n', crystal.units.length);

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
    'C_deg3', 1.334, ...
    'C_deg4', 1.750, ...
    'N', 1.073, ...
    'O', 0.837);

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

lat0 = geom.get_lattice(sys0);
fprintf('  lattice volume = %.8e bohr^3\n', lat0.volume);
fprintf('  ||H*G - 2piI|| = %.3e\n', norm(lat0.H * lat0.G - 2*pi*eye(3), 'fro'));

Lmin0 = geom.shortest_lattice_translation(lat0.H);
fprintf('  Lmin           = %.8f bohr\n', Lmin0);
fprintf('  Lmin/2         = %.8f bohr\n', 0.5 * Lmin0);

if ~(cfg.field.ewald.rcut < 0.5 * Lmin0)
    error(['Periodic Ewald real-space cache requires rcut < Lmin/2.\n' ...
           '  rcut   = %.8f bohr\n' ...
           '  Lmin/2 = %.8f bohr\n' ...
           'Increase the supercell or reduce cfg.field.ewald.rcut.'], ...
        cfg.field.ewald.rcut, 0.5 * Lmin0);
end

%% 4. Select centered neighbor pair

fprintf('\n[4] Selecting centered %s shell %d pair...\n', cfg.relation, cfg.shell);

tSelect = tic;

selection = builder.select_centered_neighbor_pair(sys0, ...
    'Relation', cfg.relation, ...
    'Shell', cfg.shell, ...
    'StackAxis', cfg.stackAxis, ...
    'Direction', cfg.direction, ...
    'Verbose', cfg.verbose);

selectTime = toc(tSelect);

fprintf('  pair selection time = %.6f s\n', selectTime);

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
    error('Selected charged pair must be neutral for periodic Ewald fixed-charge field.');
end

io.assert_atomic_units(sys);

%% 6. Extract periodic polarization system

fprintf('\n[6] Extracting periodic polarization system...\n');

polsys = builder.extract_polarization_system(sys, struct('mode', 'periodic'));

io.assert_atomic_units(polsys);

fprintf('  polsys.n_sites      = %d\n', polsys.n_sites);
fprintf('  polsys.is_periodic  = %d\n', polsys.is_periodic);
fprintf('  polarizable sites   = %d\n', nnz(polsys.site_is_polarizable));
fprintf('  charged sites       = %d\n', nnz(abs(polsys.site_charge) > 0));

lat = geom.get_lattice(polsys);
Lmin = geom.shortest_lattice_translation(lat.H);

fprintf('  lattice volume      = %.8e bohr^3\n', lat.volume);
fprintf('  ||H*G - 2piI||      = %.3e\n', norm(lat.H * lat.G - 2*pi*eye(3), 'fro'));
fprintf('  Lmin                = %.8f bohr\n', Lmin);
fprintf('  Lmin/2              = %.8f bohr\n', 0.5 * Lmin);

if ~(cfg.field.ewald.rcut < 0.5 * Lmin)
    error(['Periodic Ewald real-space cache requires rcut < Lmin/2 after extraction.\n' ...
           '  rcut   = %.8f bohr\n' ...
           '  Lmin/2 = %.8f bohr'], ...
        cfg.field.ewald.rcut, 0.5 * Lmin);
end

%% 7. Compute periodic external field

fprintf('\n[7] Computing periodic Ewald external field...\n');

fieldParams = struct();
fieldParams.use_thole = cfg.field.use_thole_damping;

fieldParams.field = struct();
fieldParams.field.mode = cfg.field.mode;
fieldParams.field.exclude_self = cfg.field.exclude_self;
fieldParams.field.use_thole_damping = cfg.field.use_thole_damping;
fieldParams.field.target_mask = logical(polsys.site_is_polarizable(:));
fieldParams.field.source_mask = abs(polsys.site_charge(:)) > 0;

fieldParams.field.real_only = cfg.field.real_only;
fieldParams.field.kspace_mode = cfg.field.kspace_mode;
fieldParams.field.k_block_size = cfg.field.k_block_size;
fieldParams.field.kspace_memory_limit_gb = cfg.field.kspace_memory_limit_gb;
fieldParams.field.verbose = cfg.field.verbose;

fieldParams.field.ewald = cfg.field.ewald;

tField = tic;
Eext = calc.compute_external_field(polsys, fieldParams);
fieldTime = toc(tField);

% Direct call for diagnostics/parts.
fieldDirect = fieldParams.field;
fieldDirect = rmfield(fieldDirect, 'mode');

[~, fieldParts] = thole.induced_field_from_charges_periodic(polsys, fieldDirect);

fprintf('  Eext computed in %.6f s\n', fieldTime);
fprintf('  ||Eext||_F               = %.12e\n', norm(Eext, 'fro'));
fprintf('  ||Eext polarizable||_F   = %.12e\n', norm(Eext(polsys.site_is_polarizable, :), 'fro'));
fprintf('  ||Ereal||_F              = %.12e\n', norm(fieldParts.real, 'fro'));
fprintf('  ||Erecip||_F             = %.12e\n', norm(fieldParts.recip, 'fro'));
fprintf('  ||Esurf||_F              = %.12e\n', norm(fieldParts.surf, 'fro'));
fprintf('  field nK                 = %d\n', fieldParts.nK);
fprintf('  field real entries       = %d\n', fieldParts.nRealEntries);
fprintf('  field storage mode       = %s\n', fieldParts.storage_mode);
fprintf('  field Thole damping      = %d\n', fieldParts.use_thole_damping);
fprintf('  selected source charge   = %+ .12e\n', fieldParts.qtot);
fprintf('  ||field B_ewald||_2          = %.12e\n', norm(fieldParts.realCache.B_ewald));
fprintf('  ||field thole_delta||_2      = %.12e\n', norm(fieldParts.realCache.thole_delta));
fprintf('  thole/ewald coeff ratio      = %.12e\n', ...
    norm(fieldParts.realCache.thole_delta) / max(norm(fieldParts.realCache.B_ewald), eps));

if abs(fieldParts.qtot) > 1e-10
    error('Periodic external field selected sources are not neutral.');
end

%% 8. Prepare SCF problem

fprintf('\n[8] Preparing SCF problem...\n');

problem = thole.prepare_scf_problem(polsys, Eext, cfg.scf);

fprintf('  nPolSites            = %d\n', problem.nPolSites);
fprintf('  active vector length = %d\n', numel(problem.Eext_pol_vec));
fprintf('  ||Eext_pol_vec||     = %.12e\n', norm(problem.Eext_pol_vec));
fprintf('  ||alpha_pol_vec||    = %.12e\n', norm(problem.alpha_pol_vec));

%% 9. Build periodic solver operator

fprintf('\n[9] Building periodic solver operator...\n');

tOp = tic;

op = thole.make_polarization_operator(polsys, problem, ...
    'Mode', cfg.operator.mode, ...
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

opTime = toc(tOp);

local_print_operator_summary(op, opTime, 'solver');

%% 10. Solve SCF

fprintf('\n[10] Solving SCF with %s...\n', cfg.solver.method);

tSolve = tic;

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

solveWallTime = toc(tSolve);

local_print_solver_summary(solveInfo, solveWallTime, mu);

%% 11. Prepare energy operator

fprintf('\n[11] Preparing energy operator...\n');

% Energy can use either dense op.Tpol or matrix-free op.apply(muVec).
% Therefore, do not build a dense operator for large matrix-free runs.
energyOp = op;

fprintf('  using solver operator for energy.\n');
fprintf('  energy op.kind    = %s\n', energyOp.kind);
fprintf('  energy op.backend = %s\n', energyOp.backend);

%% 12. Energy

fprintf('\n[12] Computing active-space polarization energy...\n');

energy = calc.compute_total_energy_active_space(polsys, problem, mu, Eext, energyOp);

fprintf('\nEnergy breakdown:\n');
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

if abs(energy.stationary_consistency) > 1e-8
    warning('Stationary energy consistency is larger than expected.');
end

%% 12b. Same-geometry nonperiodic comparison

fprintf('\n[12b] Same-geometry nonperiodic comparison...\n');

npFieldParams = fieldParams;
npFieldParams.field.mode = 'nonperiodic';

% Remove periodic-only field controls.
periodicOnlyFieldNames = {'ewald', 'kspace_mode', 'k_block_size', ...
    'kspace_memory_limit_gb', 'real_only', 'verbose'};

for a = 1:numel(periodicOnlyFieldNames)
    name = periodicOnlyFieldNames{a};
    if isfield(npFieldParams.field, name)
        npFieldParams.field = rmfield(npFieldParams.field, name);
    end
end

tNpField = tic;
Eext_np = calc.compute_external_field(polsys, npFieldParams);
npFieldTime = toc(tNpField);

fprintf('  nonperiodic Eext time       = %.6f s\n', npFieldTime);
fprintf('  ||Eext periodic||_F         = %.12e\n', norm(Eext, 'fro'));
fprintf('  ||Eext nonperiodic||_F      = %.12e\n', norm(Eext_np, 'fro'));
fprintf('  ||Eext periodic-nonper||_F  = %.12e\n', norm(Eext - Eext_np, 'fro'));
fprintf('  rel Eext difference         = %.12e\n', ...
    norm(Eext - Eext_np, 'fro') / max(norm(Eext, 'fro'), eps));

problem_np = thole.prepare_scf_problem(polsys, Eext_np, cfg.scf);

tNpOp = tic;
op_np = thole.make_polarization_operator(polsys, problem_np, ...
    'Mode', 'nonperiodic', ...
    'Solver', cfg.solver.method, ...
    'Backend', 'auto', ...
    'UseThole', cfg.operator.use_thole, ...
    'Softening', cfg.operator.softening, ...
    'Rcut', cfg.operator.rcut, ...
    'UseMex', cfg.operator.use_mex, ...
    'Profile', false, ...
    'Verbose', false);
npOpTime = toc(tNpOp);

fprintf('  nonperiodic op build time   = %.6f s\n', npOpTime);
fprintf('  nonperiodic op.kind         = %s\n', op_np.kind);
fprintf('  nonperiodic op.backend      = %s\n', op_np.backend);

tNpSolve = tic;

switch lower(cfg.solver.method)
    case 'direct'
        npSolveOpts = struct();
        npSolveOpts.compute_residual = cfg.solver.compute_residual;
        [mu_np, info_np] = thole.solve_scf_direct(problem_np, op_np, npSolveOpts);

    case 'jacobi'
        npSolveOpts = struct();
        npSolveOpts.tol = cfg.scf.tol;
        npSolveOpts.max_iter = cfg.scf.maxIter;
        npSolveOpts.mixing = cfg.solver.jacobi_mixing;
        npSolveOpts.stop_metric = cfg.solver.stop_metric;
        npSolveOpts.verbose = cfg.scf.verbose;
        [mu_np, info_np] = thole.solve_scf_jacobi(problem_np, op_np, npSolveOpts);

    case 'gmres'
        npSolveOpts = struct();
        npSolveOpts.tol = cfg.scf.tol;
        npSolveOpts.max_iter = cfg.scf.maxIter;
        npSolveOpts.restart = cfg.solver.gmres_restart;
        npSolveOpts.verbose = cfg.scf.verbose;
        [mu_np, info_np] = thole.solve_scf_gmres(problem_np, op_np, npSolveOpts);

    case 'sor'
        npSolveOpts = struct();
        npSolveOpts.tol = cfg.scf.tol;
        npSolveOpts.max_iter = cfg.scf.maxIter;
        npSolveOpts.omega = cfg.solver.sor_omega;
        npSolveOpts.stop_metric = cfg.solver.stop_metric;
        npSolveOpts.residual_every = cfg.solver.sor_residual_every;
        npSolveOpts.verbose = cfg.scf.verbose;
        [mu_np, info_np] = thole.solve_scf_sor(problem_np, op_np, npSolveOpts);

    otherwise
        error('Unsupported cfg.solver.method "%s".', cfg.solver.method);
end

npSolveTime = toc(tNpSolve);

energy_np = calc.compute_total_energy_active_space( ...
    polsys, problem_np, mu_np, Eext_np, op_np);

fprintf('  nonperiodic solve time      = %.6f s\n', npSolveTime);

if isfield(info_np, 'converged')
    fprintf('  nonperiodic converged       = %d\n', info_np.converged);
end

if isfield(info_np, 'relres')
    fprintf('  nonperiodic relres          = %.12e\n', info_np.relres);
end

fprintf('  ||mu periodic||_F           = %.12e\n', norm(mu, 'fro'));
fprintf('  ||mu nonperiodic||_F        = %.12e\n', norm(mu_np, 'fro'));
fprintf('  ||mu periodic-nonper||_F    = %.12e\n', norm(mu - mu_np, 'fro'));
fprintf('  rel mu difference           = %.12e\n', ...
    norm(mu - mu_np, 'fro') / max(norm(mu, 'fro'), eps));

fprintf('  periodic total              = %+ .12e Ha (%+ .8f eV)\n', ...
    energy.total, energy.total * HARTREE_TO_EV);
fprintf('  nonperiodic total           = %+ .12e Ha (%+ .8f eV)\n', ...
    energy_np.total, energy_np.total * HARTREE_TO_EV);
fprintf('  periodic - nonperiodic      = %+ .12e Ha (%+ .8f eV)\n', ...
    energy.total - energy_np.total, ...
    (energy.total - energy_np.total) * HARTREE_TO_EV);

fprintf('\n  component comparison, eV:\n');
fprintf('    %-24s %14s %14s %14s\n', 'component', 'periodic', 'nonperiodic', 'diff');

names = {'polarization_self', 'external_charge_dipole', 'dipole_dipole', 'total'};
for ii = 1:numel(names)
    name = names{ii};
    Ep = energy.(name) * HARTREE_TO_EV;
    En = energy_np.(name) * HARTREE_TO_EV;
    fprintf('    %-24s %+14.8f %+14.8f %+14.8f\n', ...
        name, Ep, En, Ep - En);
end

%% 13. Visualize selected pair and induced dipoles

fprintf('\n[13] Plotting selected pair and induced dipoles...\n');

fig = figure('Name', sprintf('Periodic polarization: %s', cfg.solver.method));
ax = axes(fig);

plotTitle = sprintf('Periodic %s: E_{pol} = %.4f eV', ...
    cfg.solver.method, energy.total * HARTREE_TO_EV);

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
fprintf('  ref molecule ID        = %d\n', plotOut.ref_mol_id);
fprintf('  neighbor molecule ID   = %d\n', plotOut.neighbor_mol_id);
fprintf('  pair distance / bohr   = %.6f\n', plotOut.pair_distance);
fprintf('  pair midpoint / bohr:\n');
disp(plotOut.pair_midpoint);

%% 14. Run summary

fprintf('\nRun summary:\n');
fprintf('  file                  = %s\n', cfg.filename);
fprintf('  supercell             = [%d %d %d]\n', cfg.supercellSize);
fprintf('  relation/shell        = %s / %d\n', cfg.relation, cfg.shell);
fprintf('  ref/nbr               = %d / %d\n', selection.reference_mol_id, selection.neighbor_mol_id);
fprintf('  pair distance         = %.8f bohr\n', pairDistance0);
fprintf('  charges               = [%+.3f %+.3f]\n', cfg.pairCharges);
fprintf('  field mode            = %s\n', cfg.field.mode);
fprintf('  field Thole damped    = %d\n', cfg.field.use_thole_damping);
fprintf('  field nK              = %d\n', fieldParts.nK);
fprintf('  field storage mode    = %s\n', fieldParts.storage_mode);
fprintf('  operator mode         = %s\n', op.mode);
fprintf('  operator Thole        = %d\n', cfg.operator.use_thole);
fprintf('  solver                = %s\n', cfg.solver.method);
fprintf('  solver op.kind        = %s\n', op.kind);
fprintf('  solver op.backend     = %s\n', op.backend);
fprintf('  alpha                 = %.6g\n', cfg.operator.alpha);
fprintf('  rcut                  = %.6g bohr\n', cfg.operator.rcut);
fprintf('  kcut                  = %.6g bohr^-1\n', cfg.operator.kcut);
fprintf('  boundary              = %s\n', cfg.operator.boundary);

if isfield(solveInfo, 'rowcache_fast_path_type')
    fprintf('  rowcache fast path    = %s\n', solveInfo.rowcache_fast_path_type);
end

fprintf('  Epol                  = %+ .8f eV\n', energy.total * HARTREE_TO_EV);

fprintf('\nWorkflow completed successfully.\n');

%% Local helpers

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
fprintf('  solve wall time      = %.6f s\n', wallTime);

if isfield(info, 'solve_time')
    fprintf('  solver internal time = %.6f s\n', info.solve_time);
end

if isfield(info, 'setup_time')
    fprintf('  setup time           = %.6f s\n', info.setup_time);
end

if isfield(info, 'residual_time')
    fprintf('  residual time        = %.6f s\n', info.residual_time);
end

if isfield(info, 'final_residual_time')
    fprintf('  final residual time  = %.6f s\n', info.final_residual_time);
end

if isfield(info, 'method')
    fprintf('  method               = %s\n', info.method);
end

if isfield(info, 'tol') && ~isempty(info.tol)
    fprintf('  tolerance            = %.12e\n', info.tol);
end

if isfield(info, 'stop_metric')
    fprintf('  stop metric          = %s\n', info.stop_metric);
end

if isfield(info, 'stop_value') && ~isempty(info.stop_value) && isfinite(info.stop_value)
    fprintf('  stop value           = %.12e\n', info.stop_value);
end

if isfield(info, 'residual_every')
    if isinf(info.residual_every)
        fprintf('  residual every       = Inf/final only\n');
    else
        fprintf('  residual every       = %d\n', info.residual_every);
    end
end

if isfield(info, 'iterations')
    fprintf('  iterations           = %d\n', info.iterations);
end

if isfield(info, 'max_iter')
    fprintf('  max iterations       = %d\n', info.max_iter);
end

if isfield(info, 'max_iter_requested') && info.max_iter_requested ~= info.max_iter
    fprintf('  max iter requested   = %d\n', info.max_iter_requested);
end

if isfield(info, 'mixing')
    fprintf('  mixing               = %.6f\n', info.mixing);
end

if isfield(info, 'omega')
    fprintf('  omega                = %.6f\n', info.omega);
end

if isfield(info, 'restart')
    if isempty(info.restart)
        fprintf('  gmres restart        = [] unrestarted\n');
    else
        fprintf('  gmres restart        = %d\n', info.restart);
    end
end

if isfield(info, 'iter')
    fprintf('  gmres iter           = [%s]\n', num2str(info.iter));
end

if isfield(info, 'flag')
    fprintf('  gmres flag           = %d\n', info.flag);
end

if isfield(info, 'gmres_relres') && ~isempty(info.gmres_relres) && isfinite(info.gmres_relres)
    fprintf('  gmres relres         = %.12e\n', info.gmres_relres);
end

if isfield(info, 'converged') && ~isempty(info.converged)
    fprintf('  converged            = %d\n', logical(info.converged));
end

if isfield(info, 'relres') && ~isempty(info.relres) && isfinite(info.relres)
    fprintf('  SCF relres           = %.12e\n', info.relres);
else
    fprintf('  SCF relres           = skipped/NaN\n');
end

if isfield(info, 'max_dmu') && ~isempty(info.max_dmu) && isfinite(info.max_dmu)
    fprintf('  max dmu              = %.12e\n', info.max_dmu);
end

if isfield(info, 'delta') && ~isempty(info.delta) && isfinite(info.delta)
    fprintf('  relative delta       = %.12e\n', info.delta);
end

if isfield(info, 'operator_backend')
    fprintf('  solver op backend    = %s\n', info.operator_backend);
end

if isfield(info, 'used_rowcache_fast_path')
    fprintf('  rowcache fast path   = %d\n', info.used_rowcache_fast_path);
end

if isfield(info, 'used_periodic_fast_path')
    fprintf('  periodic fast path   = %d\n', info.used_periodic_fast_path);
end

if isfield(info, 'rowcache_fast_path_type')
    fprintf('  rowcache fast type   = %s\n', info.rowcache_fast_path_type);
end

fprintf('  ||mu||_F             = %.12e\n', norm(mu, 'fro'));
fprintf('  max |mu_i|           = %.12e\n', max(vecnorm(mu, 2, 2)));
end

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