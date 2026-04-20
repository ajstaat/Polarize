%% periodic_charged_pair_manual_ewald_convergence_test
% Manual Ewald convergence study for the periodic charged-pair calculation.
%
% This script:
%   1) builds the same charged-pair crystal benchmark system
%   2) computes the CT-style external field from the assigned charges
%   3) runs periodic direct solves for a manual sweep of
%        alpha, rcut, kcut
%   4) reports total energy and decomposition
%
% Goal:
%   Find a plateau in the periodic direct energy that is insensitive to
%   the Ewald splitting parameters, rather than relying on
%   choose_ewald_params(...).

clear; clc; close all;

fprintf('=== periodic charged-pair manual Ewald convergence test ===\n');

%% ------------------------------------------------------------------------
% User controls

cfg = struct();

% Input structure
cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename   = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

% Crystal/system construction
cfg.supercellSize = [2 5 1];
cfg.bondScale     = 1.20;

% Optional nonperiodic cutoff for building the source field
cfg.rcut_bohr = 15.0;

% Charged-pair setup
cfg.pairCharges = [+1 -1];

% Field / SCF controls
cfg.softening = 0.0;
cfg.verbose = true;

% Fixed periodic setup
cfg.scale = 1.0;
cfg.boundary = 'tinfoil';

% -------------------------------------------------------------------------
% Manual Ewald sweep
%
% Start modestly. If this is too expensive, shorten these lists.
% If you start to see stabilization, refine around that region.

cfg.alpha_list = [0.30, 0.40, 0.50];
cfg.rcut_list  = [10.0, 12.0];
cfg.kcut_list  = [3.0, 3.5];

tol = struct();
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

sys0 = builder.make_crystal_system(crystal, model, buildOpts);

fprintf('\nBuilt working system:\n');
fprintf('  n_unit_sites = %d\n', sys0.n_unit_sites);
fprintf('  n_cells      = %d\n', sys0.n_cells);
fprintf('  n_sites      = %d\n', sys0.n_sites);

%% ------------------------------------------------------------------------
% Automatically choose charged pair

[refID, refSummary] = builder.choose_center_reference_molecule(sys0, ...
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

iRef = find(refSummary.candidate_mol_ids == refID, 1, 'first');
if ~isempty(iRef)
    fprintf('  ref distance to center = %.4f\n', refSummary.candidate_distance(iRef));
end

%% ------------------------------------------------------------------------
% Apply charges and disable polarizability on charged molecules

sys = builder.apply_molecule_charges(sys0, [refID nbrID], ...
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
% Canonical polarization system

params_extract = struct();
params_extract.ewald = struct();
params_extract.ewald.mode = 'periodic';

polsys0 = builder.extract_polarization_system(sys, params_extract);
io.assert_atomic_units(polsys0);

fprintf('\nCanonical polarization system:\n');
fprintf('  n_sites             = %d\n', polsys0.n_sites);
fprintf('  n_polarizable_sites = %d\n', nnz(polsys0.site_is_polarizable));

%% ------------------------------------------------------------------------
% External field from charges

fieldParams = struct();
fieldParams.field = struct();
fieldParams.field.include_external_charges = true;
fieldParams.field.exclude_self = true;
fieldParams.field.softening = cfg.softening;
if ~isempty(cfg.rcut_bohr)
    fieldParams.field.rcut = cfg.rcut_bohr;
end

fieldCacheOpts = struct();
fieldCacheOpts.site_mask = true(polsys0.n_sites, 1);

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

fieldParams.field.geom_cache = geom.build_nonperiodic_pair_cache(polsys0, fieldCacheOpts);

tField = tic;
Eext0 = calc.compute_external_field(polsys0, fieldParams);
time_external_field = toc(tField);

if isempty(Eext0) || ~isequal(size(Eext0), [polsys0.n_sites, 3])
    error('External field was not built correctly.');
end
if any(~isfinite(Eext0(:)))
    error('External field contains non-finite values.');
end

fprintf('\nExternal field:\n');
fprintf('  build time = %.6f s\n', time_external_field);

%% ------------------------------------------------------------------------
% Shared direct SCF controls

scfBase = struct();
scfBase.softening = cfg.softening;
scfBase.tol = 1e-10;
scfBase.maxIter = 500;
scfBase.verbose = cfg.verbose;
scfBase.printEvery = 25;
scfBase.residualEvery = 25;
scfBase.stopMetric = 'max_dmu';
scfBase.use_thole = false;
if ~isempty(cfg.rcut_bohr)
    scfBase.rcut = cfg.rcut_bohr;
end

%% ------------------------------------------------------------------------
% Fixed geometry / scale-1 periodic system

H0 = local_get_direct_lattice(polsys0);
polsys = polsys0;
polsys.super_lattice = H0 * cfg.scale;

center0 = 0.5 * sum(H0, 2);
center_per = 0.5 * sum(polsys.super_lattice, 2);
shift = (center_per - center0).';
polsys.site_pos = polsys0.site_pos + shift;

Eext = Eext0;
problem = thole.prepare_scf_problem(polsys, Eext, scfBase);

%% ------------------------------------------------------------------------
% Manual 3-parameter sweep

fprintf('\n--- manual Ewald sweep ---\n');

nAlpha = numel(cfg.alpha_list);
nRcut  = numel(cfg.rcut_list);
nKcut  = numel(cfg.kcut_list);

nCases = nAlpha * nRcut * nKcut;

results = struct();
results.alpha = nan(nCases, 1);
results.rcut = nan(nCases, 1);
results.kcut = nan(nCases, 1);
results.relres = nan(nCases, 1);
results.energy = nan(nCases, 1);
results.Uint = nan(nCases, 1);
results.Ureal = nan(nCases, 1);
results.Urecip = nan(nCases, 1);
results.Uself = nan(nCases, 1);
results.Usurf = nan(nCases, 1);
results.assembly_time = nan(nCases, 1);
results.solve_time = nan(nCases, 1);
results.nRealInteractions = nan(nCases, 1);
results.num_kvec = nan(nCases, 1);

caseIdx = 0;

for ia = 1:nAlpha
    alpha = cfg.alpha_list(ia);

    for ir = 1:nRcut
        rcut = cfg.rcut_list(ir);

        for ik = 1:nKcut
            kcut = cfg.kcut_list(ik);

            caseIdx = caseIdx + 1;

            fprintf('\n>>> case %d / %d : alpha = %.4f | rcut = %.4f | kcut = %.4f\n', ...
                caseIdx, nCases, alpha, rcut, kcut);

            ewaldParams = struct();
            ewaldParams.alpha = alpha;
            ewaldParams.rcut = rcut;
            ewaldParams.kcut = kcut;
            ewaldParams.boundary = cfg.boundary;

            tAsm = tic;
            [Tper, parts_per, opinfo_per] = ewald.assemble_periodic_interaction_matrix( ...
                polsys, problem, ewaldParams, scfBase);
            results.assembly_time(caseIdx) = toc(tAsm);

            tSolve = tic;
            [mu_per, direct_per] = thole.solve_scf_direct(problem, Tper);
            results.solve_time(caseIdx) = toc(tSolve);

            relres_per = thole.compute_active_space_relres(problem, Tper, mu_per);
            energy_per = calc.compute_total_energy_active_space(polsys, problem, mu_per, Eext, Tper);

            mu_per_pol = mu_per(problem.activeSites, :);
            mu_per_vec = util.stack_xyz(mu_per_pol);

            Uint_per = 0.5 * (mu_per_vec.' * (Tper * mu_per_vec));
            Ureal  = 0.5 * (mu_per_vec.' * (parts_per.real  * mu_per_vec));
            Urecip = 0.5 * (mu_per_vec.' * (parts_per.recip * mu_per_vec));
            Uself  = 0.5 * (mu_per_vec.' * (parts_per.self  * mu_per_vec));
            Usurf  = 0.5 * (mu_per_vec.' * (parts_per.surf  * mu_per_vec));

            results.alpha(caseIdx) = alpha;
            results.rcut(caseIdx) = rcut;
            results.kcut(caseIdx) = kcut;
            results.relres(caseIdx) = relres_per;
            results.energy(caseIdx) = energy_per.total;
            results.Uint(caseIdx) = Uint_per;
            results.Ureal(caseIdx) = Ureal;
            results.Urecip(caseIdx) = Urecip;
            results.Uself(caseIdx) = Uself;
            results.Usurf(caseIdx) = Usurf;
            results.nRealInteractions(caseIdx) = opinfo_per.nRealInteractions;
            results.num_kvec(caseIdx) = opinfo_per.num_kvec;

            fprintf('  periodic relres           = %.3e\n', relres_per);
            fprintf('  periodic total energy     = %+0.12f Ha\n', energy_per.total);
            fprintf('  periodic 0.5*mu''Tmu       = %+0.12f Ha\n', Uint_per);
            fprintf('    real-space contribution = %+0.12f Ha\n', Ureal);
            fprintf('    reciprocal contribution = %+0.12f Ha\n', Urecip);
            fprintf('    self contribution       = %+0.12f Ha\n', Uself);
            fprintf('    surface contribution    = %+0.12f Ha\n', Usurf);
            fprintf('    sum(parts)              = %+0.12f Ha\n', Ureal + Urecip + Uself + Usurf);
            fprintf('  nRealInteractions         = %d\n', opinfo_per.nRealInteractions);
            fprintf('  num_kvec                  = %d\n', opinfo_per.num_kvec);

            if ~isfinite(relres_per) || relres_per > tol.relres_direct
                error('Periodic direct residual too large for alpha=%.4f rcut=%.4f kcut=%.4f: %.3e', ...
                    alpha, rcut, kcut, relres_per);
            end
        end
    end
end

%% ------------------------------------------------------------------------
% Summary

fprintf('\n--- manual Ewald convergence summary ---\n');
fprintf('%8s  %8s  %8s  %14s  %14s  %14s  %10s  %10s\n', ...
    'alpha', 'rcut', 'kcut', 'Etotal(Ha)', 'Uint(Ha)', 'Urecip(Ha)', 'nReal', 'nK');

for i = 1:nCases
    fprintf('%8.4f  %8.3f  %8.3f  %14.6e  %14.6e  %14.6e  %10d  %10d\n', ...
        results.alpha(i), results.rcut(i), results.kcut(i), ...
        results.energy(i), results.Uint(i), results.Urecip(i), ...
        results.nRealInteractions(i), results.num_kvec(i));
end

% Best-effort plateau hints
fprintf('\nSimple nearest-neighbor energy deltas in scan order:\n');
for i = 2:nCases
    dE = results.energy(i) - results.energy(i-1);
    fprintf('  case %2d -> %2d : dE = %+0.12e Ha\n', i-1, i, dE);
end

fprintf('\nInterpretation:\n');
fprintf('  Look for a region where Etotal changes only weakly as rcut and kcut\n');
fprintf('  are increased for a fixed alpha, and where nearby alpha values agree.\n');
fprintf('Done.\n');

%% ========================================================================
% Local helpers
%% ========================================================================

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('local_get_direct_lattice:MissingLattice', ...
            'Need sys.super_lattice or sys.lattice.');
    end
end