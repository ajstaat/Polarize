%% periodic_operator_consistency_toy_test
% Compare assembled periodic operator pieces against matrix-free periodic applies
% on a small toy system so debugging is fast.

clear; clc; close all;

fprintf('=== periodic operator consistency toy test ===\n');

%% ------------------------------------------------------------------------
% User controls

cfg = struct();

cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename   = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

% Keep this small so the test is fast.
cfg.supercellSize = [2 5 1];
cfg.bondScale     = 1.20;

cfg.pairCharges = [+1 -1];

cfg.use_thole = true;
cfg.verbose = true;
cfg.softening = 0.0;

cfg.periodic = struct();
cfg.periodic.alpha    = 0.30;
cfg.periodic.rcut     = 12.0;
cfg.periodic.kcut     = 3.5;
cfg.periodic.boundary = 'tinfoil';

% Trial dipole choice:
%   'random'      -> random active-space dipoles
%   'direct_mu'   -> converged direct-solve dipoles
cfg.mu_test_mode = 'random';
cfg.random_seed = 1;
cfg.random_scale = 1e-3;   % a.u.

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
% Periodic params

ewaldParams = struct();
ewaldParams.alpha    = cfg.periodic.alpha;
ewaldParams.rcut     = cfg.periodic.rcut;
ewaldParams.kcut     = cfg.periodic.kcut;
ewaldParams.boundary = cfg.periodic.boundary;

scfParams = struct();
scfParams.use_thole = cfg.use_thole;
scfParams.softening = cfg.softening;
scfParams.verbose = true;
scfParams.printEvery = 250;
scfParams.residualEvery = 10;
scfParams.tol = 1e-10;
scfParams.maxIter = 500;
scfParams.omega = 0.97;
scfParams.stopMetric = 'relres';

fprintf('\n============================================================\n');
fprintf('PERIODIC TOY OPERATOR TEST\n');
fprintf('============================================================\n');
fprintf('  alpha    = %.6f\n', ewaldParams.alpha);
fprintf('  rcut     = %.6f bohr\n', ewaldParams.rcut);
fprintf('  kcut     = %.6f bohr^-1\n', ewaldParams.kcut);
fprintf('  boundary = %s\n', ewaldParams.boundary);

%% ------------------------------------------------------------------------
% Build periodic external field (needed only if we want direct_mu mode)

params_per = struct();
params_per.use_thole = cfg.use_thole;
params_per.field = struct();
params_per.field.mode = 'periodic';
params_per.field.exclude_self = true;
params_per.field.use_thole_damping = cfg.use_thole;
params_per.field.target_mask = targetMask;
params_per.field.source_mask = sourceMask;
params_per.field.ewald = ewaldParams;

fprintf('[toy] building periodic external field...\n');
Eext = calc.compute_external_field(polsys, params_per);
fprintf('[toy] ||Eext||_F = %.16e\n', norm(Eext, 'fro'));

fprintf('[toy] preparing SCF problem...\n');
problem = thole.prepare_scf_problem(polsys, Eext, scfParams);
fprintf('[toy] nPolSites = %d\n', problem.nPolSites);

%% ------------------------------------------------------------------------
% Assemble reference Tpol and component pieces

fprintf('[toy] assembling periodic operator...\n');
[Tpol, parts, opinfo] = ewald.assemble_periodic_interaction_matrix( ...
    polsys, problem, ewaldParams, scfParams); %#ok<NASGU>

fprintf('[toy] assembled operator size = %d x %d\n', size(Tpol,1), size(Tpol,2));
if isfield(opinfo, 'nRealInteractions')
    fprintf('[toy] nRealInteractions = %d\n', opinfo.nRealInteractions);
end
if isfield(opinfo, 'num_kvec')
    fprintf('[toy] num_kvec = %d\n', opinfo.num_kvec);
end

%% ------------------------------------------------------------------------
% Build matrix-free periodic caches used by iterative path

cacheParams = struct();
cacheParams.use_thole = cfg.use_thole;

fprintf('[toy] building matrix-free periodic caches...\n');
realCache_mf = geom.build_periodic_realspace_cache(polsys, problem, ewaldParams, cacheParams);
rowCache_mf  = geom.build_periodic_realspace_row_cache(polsys, problem, realCache_mf);
kCache_mf    = geom.build_periodic_kspace_cache(polsys, problem, ewaldParams);

dipoleParams = struct();
dipoleParams.use_thole = cfg.use_thole;
dipoleParams.problem = problem;
dipoleParams.target_mask = problem.polMask;
dipoleParams.source_mask = problem.polMask;
dipoleParams.realspace_cache = realCache_mf;
dipoleParams.kspace_cache = kCache_mf;

%% ------------------------------------------------------------------------
% Choose trial dipoles

activeSites = problem.activeSites(:);
nPol = problem.nPolSites;

mu_test = zeros(problem.nSites, 3);

switch lower(cfg.mu_test_mode)
    case 'random'
        rng(cfg.random_seed);
        mu_pol = cfg.random_scale * randn(nPol, 3);
        mu_test(activeSites, :) = mu_pol;

    case 'direct_mu'
        fprintf('[toy] solving direct SCF to get mu_test...\n');
        [mu_direct_ref, ~] = thole.solve_scf_direct(problem, Tpol);
        mu_test = mu_direct_ref;

    otherwise
        error('Unknown cfg.mu_test_mode = %s', cfg.mu_test_mode);
end

mu_test_pol = mu_test(activeSites, :);

fprintf('[toy] ||mu_test||_2 = %.16e\n', norm(mu_test));

%% ------------------------------------------------------------------------
% Full assembled action

mu_vec = reshape(mu_test_pol.', [], 1);
Edip_Tpol_vec = Tpol * mu_vec;
Edip_Tpol_pol = reshape(Edip_Tpol_vec, 3, []).';
Edip_Tpol = zeros(problem.nSites, 3);
Edip_Tpol(activeSites, :) = Edip_Tpol_pol;

%% ------------------------------------------------------------------------
% Full matrix-free action

Edip_mf = thole.induced_field_from_dipoles_thole_periodic( ...
    polsys, mu_test, ewaldParams, dipoleParams);

%% ------------------------------------------------------------------------
% Real-space matrix-free action only

Edip_real_mf = local_apply_periodic_realspace_only(problem, mu_test_pol, rowCache_mf);

%% ------------------------------------------------------------------------
% Reciprocal + self + surface matrix-free remainder

Edip_rest_mf = Edip_mf;
Edip_rest_mf(activeSites, :) = Edip_rest_mf(activeSites, :) - Edip_real_mf;

%% ------------------------------------------------------------------------
% Assembled component actions

Edip_real_T_vec  = parts.real  * mu_vec;
Edip_recip_T_vec = parts.recip * mu_vec;
Edip_self_T_vec  = parts.self  * mu_vec;
Edip_surf_T_vec  = parts.surf  * mu_vec;

Edip_real_T_pol  = reshape(Edip_real_T_vec,  3, []).';
Edip_recip_T_pol = reshape(Edip_recip_T_vec, 3, []).';
Edip_self_T_pol  = reshape(Edip_self_T_vec,  3, []).';
Edip_surf_T_pol  = reshape(Edip_surf_T_vec,  3, []).';

Edip_rest_T_pol = Edip_recip_T_pol + Edip_self_T_pol + Edip_surf_T_pol;

%% ------------------------------------------------------------------------
% Compare actions

diff_full_abs  = norm(Edip_mf(activeSites,:) - Edip_Tpol(activeSites,:), 'fro');
diff_full_rel  = diff_full_abs / max(norm(Edip_Tpol(activeSites,:), 'fro'), 1.0);

diff_real_abs  = norm(Edip_real_mf - Edip_real_T_pol, 'fro');
diff_real_rel  = diff_real_abs / max(norm(Edip_real_T_pol, 'fro'), 1.0);

diff_rest_abs  = norm(Edip_rest_mf(activeSites,:) - Edip_rest_T_pol, 'fro');
diff_rest_rel  = diff_rest_abs / max(norm(Edip_rest_T_pol, 'fro'), 1.0);

fprintf('\n============================================================\n');
fprintf('OPERATOR ACTION COMPARISON\n');
fprintf('============================================================\n');

fprintf('Full operator:\n');
fprintf('  ||Edip_mf - Edip_Tpol||_F      = %.16e\n', diff_full_abs);
fprintf('  relative difference            = %.16e\n', diff_full_rel);

fprintf('\nReal-space only:\n');
fprintf('  ||Edip_real_mf - Edip_real_T|| = %.16e\n', diff_real_abs);
fprintf('  relative difference            = %.16e\n', diff_real_rel);

fprintf('\nRecip+self+surf remainder:\n');
fprintf('  ||rest_mf - rest_T||_F         = %.16e\n', diff_rest_abs);
fprintf('  relative difference            = %.16e\n', diff_rest_rel);

%% ------------------------------------------------------------------------
% Compare diagonal blocks assembled vs SOR-style construction

fprintf('\n============================================================\n');
fprintf('DIAGONAL BLOCK CHECK\n');
fprintf('============================================================\n');

Ddiag_sor = local_build_periodic_sor_Ddiag(polsys, problem, ewaldParams, scfParams);

nShow = min(5, nPol);
for i = 1:nShow
    Ii = util.block3(i);
    Tii = Tpol(Ii, Ii);
    Dii = Ddiag_sor(:,:,i);

    blockDiffAbs = norm(Tii - Dii, 'fro');
    blockDiffRel = blockDiffAbs / max(norm(Tii, 'fro'), 1.0);

    fprintf('site %d (full site %d): ||Tii - Dii||_F = %.16e | rel = %.16e\n', ...
        i, activeSites(i), blockDiffAbs, blockDiffRel);
end

fprintf('\nDone.\n');

%% ========================================================================
% Local helpers
%% ========================================================================

function Edip_pol = local_apply_periodic_realspace_only(problem, mu_pol, rowCache)
    nPol = problem.nPolSites;

    row_ptr = rowCache.row_ptr;
    col_idx = rowCache.col_idx;
    dr = rowCache.dr;
    coeff_iso = rowCache.coeff_iso(:);
    coeff_dyad = rowCache.coeff_dyad(:);

    Edip_pol = zeros(nPol, 3);

    for i = 1:nPol
        idx0 = row_ptr(i);
        idx1 = row_ptr(i + 1) - 1;

        if idx1 < idx0
            continue;
        end

        idx = idx0:idx1;
        cols = col_idx(idx);
        muNbr = mu_pol(cols, :);
        rij = dr(idx, :);

        muDotR = sum(muNbr .* rij, 2);
        contrib = coeff_iso(idx) .* muNbr + coeff_dyad(idx) .* (muDotR .* rij);

        Edip_pol(i, :) = sum(contrib, 1);
    end
end

function Ddiag = local_build_periodic_sor_Ddiag(sys, problem, ewaldParams, scfParams)
    % Rebuild the same Ddiag logic used in solve_scf_iterative_periodic_sor.

    use_thole = true;
    if isfield(scfParams, 'use_thole') && ~isempty(scfParams.use_thole)
        use_thole = logical(scfParams.use_thole);
    end

    boundary = 'tinfoil';
    if isfield(ewaldParams, 'boundary') && ~isempty(ewaldParams.boundary)
        boundary = lower(char(string(ewaldParams.boundary)));
    end

    cacheParams = struct();
    cacheParams.use_thole = use_thole;

    realCache = geom.build_periodic_realspace_cache(sys, problem, ewaldParams, cacheParams);
    rowCache  = geom.build_periodic_realspace_row_cache(sys, problem, realCache);
    kCache    = geom.build_periodic_kspace_cache(sys, problem, ewaldParams);

    nPolSites = problem.nPolSites;
    activeSites = problem.activeSites(:);
    alpha_pol = problem.alpha_pol(:);

    I3 = eye(3);

    row_ptr = rowCache.row_ptr;
    col_idx = rowCache.col_idx;
    dr_all = rowCache.dr;
    coeff_iso_all = rowCache.coeff_iso(:);
    coeff_dyad_all = rowCache.coeff_dyad(:);

    % Real-space self-image diagonal blocks
    Dreal_diag = zeros(3, 3, nPolSites);

    for i = 1:nPolSites
        idx0 = row_ptr(i);
        idx1 = row_ptr(i + 1) - 1;

        if idx1 < idx0
            continue;
        end

        idx = idx0:idx1;
        selfMask = (col_idx(idx) == i);

        if any(selfMask)
            idxs = idx(selfMask);
            Dii = zeros(3,3);

            for p = idxs
                x = dr_all(p, :).';
                Dii = Dii + coeff_iso_all(p) * I3 + coeff_dyad_all(p) * (x * x.');
            end

            Dreal_diag(:,:,i) = Dii;
        end
    end

    % Reciprocal-space diagonal block
    Drecip_diag = zeros(3,3);
    if kCache.num_kvec > 0
        two_pref = 2 * kCache.pref(:);
        kvecs = kCache.kvecs;
        for m = 1:kCache.num_kvec
            k = kvecs(m, :).';
            Drecip_diag = Drecip_diag + two_pref(m) * (k * k.');
        end
    end

    alpha_ewald = ewaldParams.alpha;
    self_coeff = (4 * alpha_ewald^3 / (3 * sqrt(pi)));
    Dself = self_coeff * I3;

    surf_coeff = 0.0;
    switch boundary
        case 'tinfoil'
            surf_coeff = 0.0;
        case 'vacuum'
            H = local_get_direct_lattice(sys);
            V = abs(det(H));
            surf_coeff = 4 * pi / (3 * V);
        otherwise
            error('Unknown boundary');
    end
    Dsurf = surf_coeff * I3;

    Ddiag = zeros(3,3,nPolSites);
    for i = 1:nPolSites
        Ddiag(:,:,i) = Dreal_diag(:,:,i) + Drecip_diag + Dself + Dsurf;
        % Msolve would be: eye(3) - alpha_pol(i)*Ddiag(:,:,i)
        %#ok<NASGU>
    end
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