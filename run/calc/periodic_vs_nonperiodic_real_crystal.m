%% periodic_alpha_sweep_real_crystal
% Periodic-only polarization calculation on a real charged-pair crystal,
% sweeping Ewald alpha for a fixed [2 5 1] supercell.

clear; clc; close all;

fprintf('=== periodic alpha sweep (real crystal) ===\n');

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

% Sweep values
cfg.alpha_list = [0.20 0.25 0.30 0.35 0.40];

% Fixed periodic settings during alpha sweep
cfg.periodic = struct();
cfg.periodic.rcut  = 12.0;
cfg.periodic.kcut  = 3.5;
cfg.periodic.boundary = 'tinfoil';

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
% Shared SCF controls

scfParams = struct();
scfParams.use_thole = cfg.use_thole;
scfParams.softening = cfg.softening;
scfParams.verbose = true;
scfParams.printEvery = 250;

%% ------------------------------------------------------------------------
% Sweep

nA = numel(cfg.alpha_list);

summary = table('Size', [nA 12], ...
    'VariableTypes', { ...
        'double','double','double','double','double','double', ...
        'double','double','double','double','double','double'}, ...
    'VariableNames', { ...
        'alpha','rcut','kcut', ...
        'field_time_s','problem_time_s','operator_time_s','solve_time_s', ...
        'total_Ha','self_Ha','extchg_Ha','dd_Ha','relres'});

for ia = 1:nA
    alpha = cfg.alpha_list(ia);

    fprintf('\n============================================================\n');
    fprintf('SWEEP CASE %d / %d : alpha = %.4f\n', ia, nA, alpha);
    fprintf('============================================================\n');

    ewaldParams = struct();
    ewaldParams.alpha = alpha;
    ewaldParams.rcut = cfg.periodic.rcut;
    ewaldParams.kcut = cfg.periodic.kcut;
    ewaldParams.boundary = cfg.periodic.boundary;

    params_per = struct();
    params_per.use_thole = cfg.use_thole;
    params_per.field = struct();
    params_per.field.mode = 'periodic';
    params_per.field.exclude_self = true;
    params_per.field.use_thole_damping = cfg.use_thole;
    params_per.field.target_mask = targetMask;
    params_per.field.source_mask = sourceMask;
    params_per.field.ewald = ewaldParams;

    fprintf('[periodic] building external field...\n');
    tField = tic;
    Eext = calc.compute_external_field(polsys, params_per);
    timeField = toc(tField);
    fprintf('[periodic] external field done in %.6f s | ||Eext||_F = %.16e\n', ...
        timeField, norm(Eext, 'fro'));

    fprintf('[periodic] preparing SCF problem...\n');
    tProblem = tic;
    problem = thole.prepare_scf_problem(polsys, Eext, scfParams);
    timeProblem = toc(tProblem);
    fprintf('[periodic] problem prep done in %.6f s | nPolSites = %d\n', ...
        timeProblem, problem.nPolSites);

    fprintf('[periodic] assembling periodic operator...\n');
    tAsm = tic;
    [Tpol, Tparts, opinfo] = ewald.assemble_periodic_interaction_matrix( ...
        polsys, problem, ewaldParams, scfParams); %#ok<NASGU,ASGLU>
    timeAsm = toc(tAsm);
    fprintf('[periodic] operator done in %.6f s\n', timeAsm);
    if isfield(opinfo, 'nRealInteractions')
        fprintf('  nRealInteractions = %d\n', opinfo.nRealInteractions);
    end
    if isfield(opinfo, 'num_kvec')
        fprintf('  num_kvec          = %d\n', opinfo.num_kvec);
    end

    fprintf('[periodic] solving direct SCF...\n');
    tSolve = tic;
    [mu, direct] = thole.solve_scf_direct(problem, Tpol); %#ok<NASGU>
    timeSolve = toc(tSolve);
    fprintf('[periodic] direct solve done in %.6f s | ||mu||_2 = %.16e\n', ...
        timeSolve, norm(mu));

    fprintf('[periodic] computing energies...\n');
    E = calc.compute_total_energy_active_space(polsys, problem, mu, Eext, Tpol);
    relres = thole.compute_active_space_relres(problem, Tpol, mu);
    fprintf('[periodic] energy postprocessing done | relres = %.16e\n', relres);

    fprintf('Energy summary:\n');
    fprintf('  total                  = %+0.12e Ha\n', E.total);
    if isfield(E, 'polarization_self')
        fprintf('  polarization_self      = %+0.12e Ha\n', E.polarization_self);
    end
    if isfield(E, 'external_charge_dipole')
        fprintf('  external_charge_dipole = %+0.12e Ha\n', E.external_charge_dipole);
    end
    if isfield(E, 'dipole_dipole')
        fprintf('  dipole_dipole          = %+0.12e Ha\n', E.dipole_dipole);
    end

    summary.alpha(ia) = alpha;
    summary.rcut(ia) = ewaldParams.rcut;
    summary.kcut(ia) = ewaldParams.kcut;
    summary.field_time_s(ia) = timeField;
    summary.problem_time_s(ia) = timeProblem;
    summary.operator_time_s(ia) = timeAsm;
    summary.solve_time_s(ia) = timeSolve;
    summary.total_Ha(ia) = E.total;
    summary.self_Ha(ia) = local_get_or_nan(E, 'polarization_self');
    summary.extchg_Ha(ia) = local_get_or_nan(E, 'external_charge_dipole');
    summary.dd_Ha(ia) = local_get_or_nan(E, 'dipole_dipole');
    summary.relres(ia) = relres;
end

%% ------------------------------------------------------------------------
% Final summary

fprintf('\n============================================================\n');
fprintf('FINAL SUMMARY TABLE\n');
fprintf('============================================================\n');
disp(summary);

fprintf('\nAlpha sweep totals (Ha):\n');
for ia = 1:nA
    fprintf('  alpha = %.4f | total = %+0.12e Ha\n', ...
        summary.alpha(ia), summary.total_Ha(ia));
end

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

function val = local_get_or_nan(s, field)
    if isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = NaN;
    end
end