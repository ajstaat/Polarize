%% run_vasp_periodic_ewald_parameter_sweep
%
% Periodic Ewald parameter sweep for the VASP polarization workflow.
%
% Purpose:
%   Test convergence / alpha-independence of the periodic Ewald fixed-charge
%   field and induced-dipole operator.
%
% This workflow:
%   1. Imports a VASP structure.
%   2. Builds a [2 5 1] periodic supercell.
%   3. Selects a centered same-stack charged molecular pair.
%   4. Applies +1/-1 charges.
%   5. Extracts a periodic polarization system.
%   6. Sweeps alpha, rcut, kcut.
%   7. For each parameter set:
%        - computes periodic Ewald Eext
%        - builds periodic Ewald operator
%        - solves SCF
%        - computes active-space polarization energy
%        - prints diagnostics
%
% Notes:
%   Polarize uses direct lattice vectors as ROWS:
%
%       cart = frac * H
%
%   Public workflow code should not transpose the lattice manually.

clear;
clc;
close all;

fprintf('\n============================================================\n');
fprintf('Periodic Ewald parameter sweep workflow\n');
fprintf('============================================================\n');

%% ------------------------------------------------------------------------
% User controls
% -------------------------------------------------------------------------

cfg = struct();

cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

% Fixed for this requested sweep.
cfg.supercellSize = [2 5 1];

cfg.bondScale = 1.20;

cfg.relation = 'same_stack';
cfg.shell = 1;
cfg.stackAxis = 'b';
cfg.direction = 'either';

cfg.pairCharges = [+1 -1];

cfg.verbose = true;

% Sweep controls.
%
% A compact default grid that checks:
%   - alpha dependence at fixed rcut/kcut
%   - kcut convergence
%   - rcut sensitivity
%
% Expand as needed after the first pass.
cfg.sweep = struct();

cfg.sweep.alphaList = [0.20 0.25 0.30 0.35 0.40];
cfg.sweep.rcutList  = [14.0 16.0 18.0 20.0];
cfg.sweep.kcutList  = [0.75 1.00 1.25 1.50 2.00 2.50];

% If true, run full Cartesian product:
%
%   numel(alphaList) * numel(rcutList) * numel(kcutList)
%
% If false, run three smaller line sweeps around the reference point:
%   1. alpha sweep at reference rcut/kcut
%   2. rcut sweep at reference alpha/kcut
%   3. kcut sweep at reference alpha/rcut
cfg.sweep.fullCartesian = false;

cfg.sweep.refAlpha = 0.30;
cfg.sweep.refRcut  = 18.0;
cfg.sweep.refKcut  = 1.25;

% Solver controls.
cfg.solver = struct();
cfg.solver.method = 'sor';          % 'direct' | 'jacobi' | 'gmres' | 'sor'
cfg.solver.stop_metric = 'max_dmu'; % 'max_dmu' or 'relres'
cfg.solver.sor_omega = 0.97;
cfg.solver.sor_residual_every = 25;
cfg.solver.jacobi_mixing = 0.6;
cfg.solver.gmres_restart = [];

% SCF controls.
cfg.scf = struct();
cfg.scf.tol = 1e-8;
cfg.scf.maxIter = 500;
cfg.scf.mixing = 0.6;
cfg.scf.omega = cfg.solver.sor_omega;
cfg.scf.verbose = false;

% Field/operator controls.
cfg.field = struct();
cfg.field.mode = 'periodic';
cfg.field.exclude_self = true;
cfg.field.use_thole_damping = true;
cfg.field.kspace_mode = 'auto';
cfg.field.k_block_size = 2048;
cfg.field.kspace_memory_limit_gb = 8;
cfg.field.real_only = false;
cfg.field.verbose = false;
cfg.field.boundary = 'tinfoil';

cfg.operator = struct();
cfg.operator.mode = 'periodic_ewald';
cfg.operator.backend = 'auto';
cfg.operator.use_thole = true;
cfg.operator.softening = 0.0;
cfg.operator.kspace_mode = cfg.field.kspace_mode;
cfg.operator.k_block_size = cfg.field.k_block_size;
cfg.operator.kspace_memory_limit_gb = cfg.field.kspace_memory_limit_gb;
cfg.operator.use_mex = true;
cfg.operator.use_mex_kspace = true;
cfg.operator.profile = false;
cfg.operator.verbose = false;

cfg.compareNonperiodicReference = true;

HARTREE_TO_EV = 27.211386245988;

%% ------------------------------------------------------------------------
% Check file
% -------------------------------------------------------------------------

if ~isfile(cfg.filename)
    error('VASP file not found:\n  %s\nEdit cfg.filename at the top of this script.', cfg.filename);
end

fprintf('\nInput structure:\n  %s\n', cfg.filename);

fprintf('\nSweep controls:\n');
fprintf('  supercell size      = [%d %d %d]\n', cfg.supercellSize);
fprintf('  full Cartesian      = %d\n', cfg.sweep.fullCartesian);
fprintf('  reference alpha     = %.6g\n', cfg.sweep.refAlpha);
fprintf('  reference rcut      = %.6g bohr\n', cfg.sweep.refRcut);
fprintf('  reference kcut      = %.6g bohr^-1\n', cfg.sweep.refKcut);
fprintf('  boundary            = %s\n', cfg.field.boundary);
fprintf('  solver              = %s\n', cfg.solver.method);
fprintf('  field Thole damping = %d\n', cfg.field.use_thole_damping);
fprintf('  op Thole damping    = %d\n', cfg.operator.use_thole);

%% ------------------------------------------------------------------------
% 1. Import crystal template
% -------------------------------------------------------------------------

fprintf('\n[1] Importing crystal template...\n');

crystal = io.import_contcar_as_crystal(cfg.filename, ...
    'BondScale', cfg.bondScale, ...
    'SortMolecules', false);

fprintf('  nSites    = %d\n', crystal.nSites);
fprintf('  nBaseMols = %d\n', crystal.nBaseMols);
fprintf('  units     = %s\n', crystal.units.length);

%% ------------------------------------------------------------------------
% 2. Model
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

sys0 = builder.make_crystal_system(crystal, model, buildOpts);

fprintf('\nBuilt system:\n');
fprintf('  supercell_size = [%d %d %d]\n', sys0.supercell_size);
fprintf('  n_sites        = %d\n', sys0.n_sites);
fprintf('  n_molecules    = %d\n', numel(sys0.molecule_table.molecule_id));
fprintf('  n_complete     = %d\n', nnz(sys0.molecule_table.is_complete_in_display));
fprintf('  length unit    = %s\n', sys0.units.length);
fprintf('  alpha unit     = %s\n', sys0.units.alpha);

lat0 = geom.get_lattice(sys0);
Lmin0 = geom.shortest_lattice_translation(lat0.H);

fprintf('  lattice volume = %.8e bohr^3\n', lat0.volume);
fprintf('  ||H*G - 2piI|| = %.3e\n', norm(lat0.H * lat0.G - 2*pi*eye(3), 'fro'));
fprintf('  Lmin           = %.8f bohr\n', Lmin0);
fprintf('  Lmin/2         = %.8f bohr\n', 0.5 * Lmin0);

%% ------------------------------------------------------------------------
% 4. Select pair
% -------------------------------------------------------------------------

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

%% ------------------------------------------------------------------------
% 5. Apply charges
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

chargedMask = ismember(sys.site_mol_id, molIDs);

fprintf('\nCharge summary:\n');
fprintf('  ref molecule %d total charge = %+ .8f e\n', ...
    selection.reference_mol_id, sum(sys.site_charge(idxRef)));
fprintf('  nbr molecule %d total charge = %+ .8f e\n', ...
    selection.neighbor_mol_id, sum(sys.site_charge(idxNbr)));
fprintf('  selected source charge       = %+ .8f e\n', ...
    sum(sys.site_charge(idxRef)) + sum(sys.site_charge(idxNbr)));
fprintf('  total system charge          = %+ .8f e\n', sum(sys.site_charge));

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

%% ------------------------------------------------------------------------
% 6. Extract periodic polarization system
% -------------------------------------------------------------------------

fprintf('\n[6] Extracting periodic polarization system...\n');

polsys = builder.extract_polarization_system(sys, struct('mode', 'periodic'));

io.assert_atomic_units(polsys);

lat = geom.get_lattice(polsys);
Lmin = geom.shortest_lattice_translation(lat.H);

fprintf('  polsys.n_sites      = %d\n', polsys.n_sites);
fprintf('  polsys.is_periodic  = %d\n', polsys.is_periodic);
fprintf('  polarizable sites   = %d\n', nnz(polsys.site_is_polarizable));
fprintf('  charged sites       = %d\n', nnz(abs(polsys.site_charge) > 0));
fprintf('  lattice volume      = %.8e bohr^3\n', lat.volume);
fprintf('  ||H*G - 2piI||      = %.3e\n', norm(lat.H * lat.G - 2*pi*eye(3), 'fro'));
fprintf('  Lmin                = %.8f bohr\n', Lmin);
fprintf('  Lmin/2              = %.8f bohr\n', 0.5 * Lmin);

sourceMask = abs(polsys.site_charge(:)) > 0;
Mq = sum(polsys.site_charge(sourceMask) .* polsys.site_pos(sourceMask, :), 1);

fprintf('\nSource pair dipole diagnostic:\n');
fprintf('  Mq / e bohr = [% .8e % .8e % .8e]\n', Mq(1), Mq(2), Mq(3));
fprintf('  |Mq|        = %.8e e bohr\n', norm(Mq));
fprintf('  |Mq|^2/V    = %.8e\n', dot(Mq, Mq) / lat.volume);

%% ------------------------------------------------------------------------
% 7. Build sweep table
% -------------------------------------------------------------------------

fprintf('\n[7] Building sweep table...\n');

sweepRows = local_build_sweep_rows(cfg);

% Remove impossible rcut values.
safeMask = sweepRows.rcut < 0.5 * Lmin;
if any(~safeMask)
    fprintf('  dropping %d rows with rcut >= Lmin/2\n', nnz(~safeMask));
    sweepRows = sweepRows(safeMask, :);
end

nRuns = height(sweepRows);

if nRuns == 0
    error('No valid sweep rows remain after rcut < Lmin/2 filtering.');
end

fprintf('  n sweep rows = %d\n', nRuns);
disp(sweepRows);

%% ------------------------------------------------------------------------
% Optional same-geometry nonperiodic reference
% -------------------------------------------------------------------------

npReference = struct();
npReference.enabled = cfg.compareNonperiodicReference;

if cfg.compareNonperiodicReference
    fprintf('\n[8] Computing same-geometry nonperiodic reference once...\n');

    npFieldParams = struct();
    npFieldParams.use_thole = cfg.field.use_thole_damping;

    npFieldParams.field = struct();
    npFieldParams.field.mode = 'nonperiodic';
    npFieldParams.field.exclude_self = cfg.field.exclude_self;
    npFieldParams.field.use_thole_damping = cfg.field.use_thole_damping;
    npFieldParams.field.target_mask = logical(polsys.site_is_polarizable(:));
    npFieldParams.field.source_mask = abs(polsys.site_charge(:)) > 0;

    tNpField = tic;
    Eext_np = calc.compute_external_field(polsys, npFieldParams);
    npFieldTime = toc(tNpField);

    problem_np = thole.prepare_scf_problem(polsys, Eext_np, cfg.scf);

    tNpOp = tic;
    op_np = thole.make_polarization_operator(polsys, problem_np, ...
        'Mode', 'nonperiodic', ...
        'Solver', cfg.solver.method, ...
        'Backend', 'auto', ...
        'UseThole', cfg.operator.use_thole, ...
        'Softening', cfg.operator.softening, ...
        'Rcut', cfg.sweep.refRcut, ...
        'UseMex', cfg.operator.use_mex, ...
        'Profile', false, ...
        'Verbose', false);
    npOpTime = toc(tNpOp);

    tNpSolve = tic;
    [mu_np, info_np] = local_solve_problem(problem_np, op_np, cfg);
    npSolveTime = toc(tNpSolve);

    energy_np = calc.compute_total_energy_active_space( ...
        polsys, problem_np, mu_np, Eext_np, op_np);

    npReference.Eext = Eext_np;
    npReference.problem = problem_np;
    npReference.op = op_np;
    npReference.mu = mu_np;
    npReference.info = info_np;
    npReference.energy = energy_np;
    npReference.fieldTime = npFieldTime;
    npReference.opTime = npOpTime;
    npReference.solveTime = npSolveTime;

    fprintf('  nonperiodic Eext time  = %.6f s\n', npFieldTime);
    fprintf('  nonperiodic op time    = %.6f s\n', npOpTime);
    fprintf('  nonperiodic solve time = %.6f s\n', npSolveTime);
    fprintf('  nonperiodic Epol       = %+ .8f eV\n', energy_np.total * HARTREE_TO_EV);
    fprintf('  nonperiodic relres     = %.12e\n', local_info_field(info_np, 'relres', NaN));
else
    fprintf('\n[8] Skipping nonperiodic reference.\n');
end

%% ------------------------------------------------------------------------
% 9. Sweep periodic Ewald parameters
% -------------------------------------------------------------------------

fprintf('\n[9] Running periodic Ewald sweep...\n');

results = table();

for irun = 1:nRuns
    alpha = sweepRows.alpha(irun);
    rcut = sweepRows.rcut(irun);
    kcut = sweepRows.kcut(irun);
    label = string(sweepRows.label(irun));

    fprintf('\n------------------------------------------------------------\n');
    fprintf('Sweep row %d / %d: %s | alpha=%.6g rcut=%.6g kcut=%.6g\n', ...
        irun, nRuns, label, alpha, rcut, kcut);
    fprintf('------------------------------------------------------------\n');

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

    fieldParams.field.ewald = struct();
    fieldParams.field.ewald.alpha = alpha;
    fieldParams.field.ewald.rcut = rcut;
    fieldParams.field.ewald.kcut = kcut;
    fieldParams.field.ewald.boundary = cfg.field.boundary;

    tField = tic;
    Eext = calc.compute_external_field(polsys, fieldParams);
    fieldTime = toc(tField);

    fieldDirect = fieldParams.field;
    fieldDirect = rmfield(fieldDirect, 'mode');
    [~, fieldParts] = thole.induced_field_from_charges_periodic(polsys, fieldDirect);

    problem = thole.prepare_scf_problem(polsys, Eext, cfg.scf);

    tOp = tic;
    op = thole.make_polarization_operator(polsys, problem, ...
        'Mode', cfg.operator.mode, ...
        'Solver', cfg.solver.method, ...
        'Backend', cfg.operator.backend, ...
        'UseThole', cfg.operator.use_thole, ...
        'Softening', cfg.operator.softening, ...
        'Rcut', rcut, ...
        'Alpha', alpha, ...
        'Kcut', kcut, ...
        'Boundary', cfg.field.boundary, ...
        'KspaceMode', cfg.operator.kspace_mode, ...
        'KBlockSize', cfg.operator.k_block_size, ...
        'KspaceMemoryLimitGB', cfg.operator.kspace_memory_limit_gb, ...
        'UseMex', cfg.operator.use_mex, ...
        'UseMexKspace', cfg.operator.use_mex_kspace, ...
        'Profile', cfg.operator.profile, ...
        'Verbose', cfg.operator.verbose);
    opTime = toc(tOp);

    tSolve = tic;
    [mu, info] = local_solve_problem(problem, op, cfg);
    solveTime = toc(tSolve);

    energy = calc.compute_total_energy_active_space(polsys, problem, mu, Eext, op);

    relres = local_info_field(info, 'relres', NaN);
    iterations = local_info_field(info, 'iterations', NaN);
    converged = local_info_field(info, 'converged', false);
    maxDmu = local_info_field(info, 'max_dmu', NaN);

    periodicFastPath = false;
    if isfield(info, 'used_periodic_fast_path')
        periodicFastPath = logical(info.used_periodic_fast_path);
    end

    opStorageMode = "";
    opEstimatedGB = NaN;
    if isfield(op, 'k_cache')
        if isfield(op.k_cache, 'storage_mode')
            opStorageMode = string(op.k_cache.storage_mode);
        end
        if isfield(op.k_cache, 'estimated_full_gb')
            opEstimatedGB = op.k_cache.estimated_full_gb;
        end
    end

    fieldEstimatedGB = NaN;
    if isfield(fieldParts, 'estimated_full_gb')
        fieldEstimatedGB = fieldParts.estimated_full_gb;
    end

    npDiffEV = NaN;
    relEextDiffNP = NaN;
    relMuDiffNP = NaN;

    if cfg.compareNonperiodicReference
        npDiffEV = (energy.total - npReference.energy.total) * HARTREE_TO_EV;
        relEextDiffNP = norm(Eext - npReference.Eext, 'fro') / max(norm(Eext, 'fro'), eps);
        relMuDiffNP = norm(mu - npReference.mu, 'fro') / max(norm(mu, 'fro'), eps);
    end

    fprintf('  field time             = %.6f s\n', fieldTime);
    fprintf('  op build time          = %.6f s\n', opTime);
    fprintf('  solve time             = %.6f s\n', solveTime);
    fprintf('  field nK               = %d\n', fieldParts.nK);
    fprintf('  operator nK            = %d\n', op.info.nK);
    fprintf('  field storage mode     = %s\n', fieldParts.storage_mode);
    fprintf('  op storage mode        = %s\n', opStorageMode);
    fprintf('  op estimated full GB   = %.3f\n', opEstimatedGB);
    fprintf('  ||Eext||_F             = %.12e\n', norm(Eext, 'fro'));
    fprintf('  ||Ereal||_F            = %.12e\n', norm(fieldParts.real, 'fro'));
    fprintf('  ||Erecip||_F           = %.12e\n', norm(fieldParts.recip, 'fro'));
    fprintf('  ||Esurf||_F            = %.12e\n', norm(fieldParts.surf, 'fro'));
    fprintf('  ||mu||_F               = %.12e\n', norm(mu, 'fro'));
    fprintf('  Epol                   = %+ .8f eV\n', energy.total * HARTREE_TO_EV);
    fprintf('  stationary consistency = %+ .12e Ha\n', energy.stationary_consistency);
    fprintf('  relres                 = %.12e\n', relres);
    fprintf('  iterations             = %g\n', iterations);
    fprintf('  converged              = %d\n', logical(converged));
    fprintf('  periodic fast path     = %d\n', periodicFastPath);

    if cfg.compareNonperiodicReference
        fprintf('  Epol - Epol_np         = %+ .8f eV\n', npDiffEV);
        fprintf('  rel Eext diff vs NP    = %.12e\n', relEextDiffNP);
        fprintf('  rel mu diff vs NP      = %.12e\n', relMuDiffNP);
    end

    row = table();
    row.label = label;
    row.alpha = alpha;
    row.rcut = rcut;
    row.kcut = kcut;

    row.nK_field = fieldParts.nK;
    row.nK_operator = op.info.nK;

    row.field_storage_mode = string(fieldParts.storage_mode);
    row.operator_storage_mode = opStorageMode;

    row.field_estimated_full_gb = fieldEstimatedGB;
    row.operator_estimated_full_gb = opEstimatedGB;

    row.field_time_s = fieldTime;
    row.operator_time_s = opTime;
    row.solve_time_s = solveTime;

    row.norm_Eext = norm(Eext, 'fro');
    row.norm_Ereal = norm(fieldParts.real, 'fro');
    row.norm_Erecip = norm(fieldParts.recip, 'fro');
    row.norm_Esurf = norm(fieldParts.surf, 'fro');
    row.norm_mu = norm(mu, 'fro');

    row.energy_self_eV = energy.polarization_self * HARTREE_TO_EV;
    row.energy_cross_eV = energy.external_charge_dipole * HARTREE_TO_EV;
    row.energy_dipdip_eV = energy.dipole_dipole * HARTREE_TO_EV;
    row.energy_total_eV = energy.total * HARTREE_TO_EV;
    row.energy_stationary_eV = energy.total_stationary * HARTREE_TO_EV;
    row.stationary_consistency_Ha = energy.stationary_consistency;
    row.energy_relres = energy.relres;

    row.solver_relres = relres;
    row.solver_iterations = iterations;
    row.solver_max_dmu = maxDmu;
    row.solver_converged = logical(converged);
    row.periodic_fast_path = periodicFastPath;

    row.diff_total_vs_np_eV = npDiffEV;
    row.rel_Eext_diff_vs_np = relEextDiffNP;
    row.rel_mu_diff_vs_np = relMuDiffNP;

    results = [results; row]; %#ok<AGROW>
end

%% ------------------------------------------------------------------------
% 10. Summary
% -------------------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('Periodic Ewald sweep summary\n');
fprintf('============================================================\n');

summaryCols = {'label', 'alpha', 'rcut', 'kcut', ...
    'nK_operator', 'operator_storage_mode', ...
    'norm_Ereal', 'norm_Erecip', 'norm_mu', ...
    'energy_total_eV', 'stationary_consistency_Ha', ...
    'solver_relres', 'solver_iterations', 'solve_time_s'};

if cfg.compareNonperiodicReference
    summaryCols = [summaryCols, {'diff_total_vs_np_eV'}];
end

disp(results(:, summaryCols));

fprintf('\nEnergy totals by row:\n');
for i = 1:height(results)
    fprintf('  %-14s alpha=%5.2f rcut=%6.2f kcut=%5.2f | Epol=%+ .8f eV | nK=%6d | relres=%.3e\n', ...
        results.label(i), results.alpha(i), results.rcut(i), results.kcut(i), ...
        results.energy_total_eV(i), results.nK_operator(i), results.solver_relres(i));
end

fprintf('\nReference nonperiodic comparison:\n');
if cfg.compareNonperiodicReference
    fprintf('  Epol_np = %+ .8f eV\n', npReference.energy.total * HARTREE_TO_EV);
else
    fprintf('  skipped\n');
end

fprintf('\nSuggested alpha-independence check:\n');
fprintf('  For rows with same rcut/kcut, energy_total_eV should be stable vs alpha.\n');
fprintf('  If not, increase rcut and/or kcut.\n');

fprintf('\nWorkflow completed successfully.\n');

%% =========================================================================
% Local helpers
% =========================================================================

function rows = local_build_sweep_rows(cfg)
if cfg.sweep.fullCartesian
    [A, R, K] = ndgrid(cfg.sweep.alphaList, cfg.sweep.rcutList, cfg.sweep.kcutList);

    alpha = A(:);
    rcut = R(:);
    kcut = K(:);
    label = strings(numel(alpha), 1);

    for i = 1:numel(alpha)
        label(i) = sprintf("grid_%03d", i);
    end

    rows = table(label, alpha, rcut, kcut);
    return;
end

label = strings(0, 1);
alpha = zeros(0, 1);
rcut = zeros(0, 1);
kcut = zeros(0, 1);

% Reference row.
label(end+1, 1) = "ref";
alpha(end+1, 1) = cfg.sweep.refAlpha;
rcut(end+1, 1) = cfg.sweep.refRcut;
kcut(end+1, 1) = cfg.sweep.refKcut;

% Alpha line at reference rcut/kcut.
for a = cfg.sweep.alphaList(:).'
    if a == cfg.sweep.refAlpha
        continue;
    end

    label(end+1, 1) = sprintf("alpha_%g", a);
    alpha(end+1, 1) = a;
    rcut(end+1, 1) = cfg.sweep.refRcut;
    kcut(end+1, 1) = cfg.sweep.refKcut;
end

% Rcut line at reference alpha/kcut.
for r = cfg.sweep.rcutList(:).'
    if r == cfg.sweep.refRcut
        continue;
    end

    label(end+1, 1) = sprintf("rcut_%g", r);
    alpha(end+1, 1) = cfg.sweep.refAlpha;
    rcut(end+1, 1) = r;
    kcut(end+1, 1) = cfg.sweep.refKcut;
end

% Kcut line at reference alpha/rcut.
for k = cfg.sweep.kcutList(:).'
    if k == cfg.sweep.refKcut
        continue;
    end

    label(end+1, 1) = sprintf("kcut_%g", k);
    alpha(end+1, 1) = cfg.sweep.refAlpha;
    rcut(end+1, 1) = cfg.sweep.refRcut;
    kcut(end+1, 1) = k;
end

rows = table(label, alpha, rcut, kcut);

% Remove accidental exact duplicates, keeping first.
[~, ia] = unique(rows(:, {'alpha','rcut','kcut'}), 'rows', 'stable');
rows = rows(ia, :);
end

function [mu, info] = local_solve_problem(problem, op, cfg)
switch lower(cfg.solver.method)
    case 'direct'
        solveOpts = struct();
        solveOpts.compute_residual = true;
        [mu, info] = thole.solve_scf_direct(problem, op, solveOpts);

    case 'jacobi'
        solveOpts = struct();
        solveOpts.tol = cfg.scf.tol;
        solveOpts.max_iter = cfg.scf.maxIter;
        solveOpts.mixing = cfg.solver.jacobi_mixing;
        solveOpts.stop_metric = cfg.solver.stop_metric;
        solveOpts.verbose = cfg.scf.verbose;
        [mu, info] = thole.solve_scf_jacobi(problem, op, solveOpts);

    case 'gmres'
        solveOpts = struct();
        solveOpts.tol = cfg.scf.tol;
        solveOpts.max_iter = cfg.scf.maxIter;
        solveOpts.restart = cfg.solver.gmres_restart;
        solveOpts.verbose = cfg.scf.verbose;
        [mu, info] = thole.solve_scf_gmres(problem, op, solveOpts);

    case 'sor'
        solveOpts = struct();
        solveOpts.tol = cfg.scf.tol;
        solveOpts.max_iter = cfg.scf.maxIter;
        solveOpts.omega = cfg.solver.sor_omega;
        solveOpts.stop_metric = cfg.solver.stop_metric;
        solveOpts.residual_every = cfg.solver.sor_residual_every;
        solveOpts.verbose = cfg.scf.verbose;
        [mu, info] = thole.solve_scf_sor(problem, op, solveOpts);

    otherwise
        error('Unsupported cfg.solver.method "%s".', cfg.solver.method);
end
end

function value = local_info_field(info, name, defaultValue)
if isstruct(info) && isfield(info, name) && ~isempty(info.(name))
    value = info.(name);
else
    value = defaultValue;
end
end