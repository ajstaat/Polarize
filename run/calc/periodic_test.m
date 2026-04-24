%% run_periodic_eext_sor_mode_comparison_with_p3m
% Controlled diagnostic:
%
%   A. nonperiodic Eext + periodic SOR/tinfoil
%   B. periodic Ewald Eext + periodic SOR/tinfoil
%   C. periodic P3M  Eext + periodic SOR/tinfoil
%   D. periodic Ewald Eext + periodic SOR/vacuum
%
% Also computes post-hoc point-dipole ghost-image corrections:
%
%   dE_mu_img  = +1/2 sum_i mu_i . E_img(Mmu)
%   dE_tot_img = +1/2 sum_i mu_i . E_img(Mq + Mmu)
%
% where
%
%   E_img(M) = E_periodic_point_dipole(M) - E_bare_central_point_dipole(M)
%
% evaluated using minimum-image displacements relative to a coherent
% charged-pair COM center.
%
% Notes:
%   - For the nonperiodic-source hybrid case, dE_mu_img is the meaningful
%     correction: the CT source images are absent, but the periodic induced
%     response repeats the induced cloud.
%
%   - For periodic-source cases, dE_tot_img is a total-cell point-dipole
%     image diagnostic.
%
%   - The P3M external field currently handles periodic fixed-charge Eext.
%     It uses:
%
%         Ereal = thole periodic real-only cache/MEX path
%         Erecip = P3M mesh reciprocal field
%         Esurf = analytic surface term, if boundary='vacuum'
%
% Assumes current fixes:
%   - thole.induced_field_from_charges_periodic supports real_only=true.
%   - p3m.compute_external_field_charges supports realspace_backend =
%     'thole_periodic_real'.
%   - thole.solve_scf_iterative_periodic_sor has the corrected periodic
%     operator/sign conventions and final residual recomputation.

clear; clc; close all;

fprintf('=== periodic Eext / SOR mode-comparison test with P3M ===\n');

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

% Ewald controls for periodic Eext and periodic SOR.
cfg.rcut = 11.5;
cfg.kcut = 3.5;
cfg.alphaList = [0.35];

% Direct Ewald Eext reciprocal controls.
cfg.field_kspace_mode = 'full';
cfg.field_k_block_size = 2048;
cfg.field_kspace_memory_limit_gb = 8;

% P3M external-field controls.
cfg.p3m = struct();
cfg.p3m.mesh_size = [32 48 32];     % speed/default
% cfg.p3m.mesh_size = [40 64 40];   % high-accuracy check
cfg.p3m.assignment_order = 4;
cfg.p3m.realspace_backend = 'thole_periodic_real';
cfg.p3m.deconvolve_assignment = true;
cfg.p3m.deconvolution_floor = 1e-6;
cfg.p3m.verbose = false;

% SOR controls.
cfg.sor = struct();
cfg.sor.use_thole = cfg.use_thole;
cfg.sor.verbose = true;
cfg.sor.printEvery = 10;
cfg.sor.residualEvery = 10;
cfg.sor.tol = 1e-6;
cfg.sor.maxIter = 500;
cfg.sor.omega = 0.55;
cfg.sor.stopMetric = 'max_dmu';

% Real-space dipole row-cache controls.
cfg.sor.use_mex_realspace = true;

% Reciprocal dipole controls inside SOR.
% Note: this is still the periodic induced-dipole reciprocal operator,
% not yet P3M for induced dipoles.
cfg.sor.kspace_mode = 'full';
cfg.sor.k_block_size = 2048;
cfg.sor.kspace_memory_limit_gb = 8;
cfg.sor.use_mex_kspace = false;

% Ghost correction controls.
% These can be the same as the case Ewald params; real-geometry debug
% showed this correction was stable across rcut/kcut/alpha after the
% minimum-image fix.
cfg.ghost.use_case_ewald_params = true;
cfg.ghost.alpha = 0.35;
cfg.ghost.rcut = 11.5;
cfg.ghost.kcut = 3.5;

% Print detailed per-site/molecule diagnostics.
cfg.printTopSites = true;
cfg.printTopMolecules = true;
cfg.nTopSites = 5;
cfg.nTopMolecules = 12;

% Cases.
cases = struct([]);

cases(1).label    = 'nonperiodic_Eext__periodic_SOR_tinfoil';
cases(1).eextMode = 'nonperiodic';
cases(1).boundary = 'tinfoil';

cases(2).label    = 'periodic_Ewald_Eext__periodic_SOR_tinfoil';
cases(2).eextMode = 'periodic_ewald';
cases(2).boundary = 'tinfoil';

cases(3).label    = 'periodic_P3M_Eext__periodic_SOR_tinfoil';
cases(3).eextMode = 'periodic_p3m';
cases(3).boundary = 'tinfoil';

cases(4).label    = 'periodic_Ewald_Eext__periodic_SOR_vacuum';
cases(4).eextMode = 'periodic_ewald';
cases(4).boundary = 'vacuum';

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

fprintf('\nChosen pair:\n');
fprintf('  ref ID = %d\n', refID);
fprintf('  nbr ID = %d\n', nbrID);

%% ------------------------------------------------------------------------
% Apply charges and extract canonical polarization system

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

polsys = builder.extract_polarization_system(sys, params_extract);
io.assert_atomic_units(polsys);

fprintf('\nCanonical polarization system:\n');
fprintf('  n_sites             = %d\n', polsys.n_sites);
fprintf('  n_polarizable_sites = %d\n', nnz(polsys.site_is_polarizable));

targetMask = logical(polsys.site_is_polarizable(:));
sourceMask = abs(polsys.site_charge(:)) > 0;

fprintf('\nMask sanity:\n');
fprintf('  n target sites        = %d\n', nnz(targetMask));
fprintf('  n source sites        = %d\n', nnz(sourceMask));
fprintf('  target/source overlap = %d\n', nnz(targetMask & sourceMask));
fprintf('  total source charge   = %+0.16e\n', sum(polsys.site_charge(sourceMask)));
fprintf('  min source alpha      = %.16e\n', min(polsys.site_alpha(sourceMask)));
fprintf('  max source alpha      = %.16e\n', max(polsys.site_alpha(sourceMask)));

if nnz(targetMask & sourceMask) ~= 0
    error('Charged source sites are still marked polarizable.');
end
if abs(sum(polsys.site_charge(sourceMask))) > 1e-10
    error('Selected periodic source charges are not neutral.');
end

%% ------------------------------------------------------------------------
% Lattice / cutoff sanity

Hrow = local_get_direct_lattice(polsys);
Hcol = Hrow.';
Lmin = geom.shortest_lattice_translation(Hcol);
rcutMaxSafe = 0.5 * Lmin - 1e-12 * max(1, Lmin);

fprintf('\nLattice/cutoff sanity:\n');
fprintf('  units.length     = %s\n', polsys.units.length);
fprintf('  volume           = %.16e bohr^3\n', abs(det(Hcol)));
fprintf('  Lmin             = %.16e bohr\n', Lmin);
fprintf('  Lmin/2           = %.16e bohr\n', 0.5 * Lmin);
fprintf('  max safe rcut    = %.16e bohr\n', rcutMaxSafe);
fprintf('  requested rcut   = %.16e bohr\n', cfg.rcut);

if ~(cfg.rcut < rcutMaxSafe)
    error('cfg.rcut violates the single-image condition for atomistic row cache.');
end

%% ------------------------------------------------------------------------
% P3M sanity print

fprintf('\nP3M external-field settings:\n');
fprintf('  mesh_size             = [%d %d %d]\n', cfg.p3m.mesh_size);
fprintf('  assignment_order      = %d\n', cfg.p3m.assignment_order);
fprintf('  realspace_backend     = %s\n', cfg.p3m.realspace_backend);
fprintf('  deconvolve_assignment = %d\n', cfg.p3m.deconvolve_assignment);
fprintf('  deconvolution_floor   = %.3e\n', cfg.p3m.deconvolution_floor);

%% ------------------------------------------------------------------------
% MEX sanity

if cfg.sor.use_mex_realspace
    fprintf('\nMEX real-space dipole row cache requested.\n');
    fprintf('  build_active_row_cache_periodic will be called with rowOpts.use_mex = true.\n');
end

if cfg.sor.use_mex_kspace
    if exist('mex_periodic_kspace_block', 'file') ~= 3
        error(['cfg.sor.use_mex_kspace=true, but mex_periodic_kspace_block is not available.\n' ...
               'Compile it first or set cfg.sor.use_mex_kspace=false.']);
    end
end

%% ------------------------------------------------------------------------
% Charged-pair ghost center

ghostCenter = local_charged_pair_com_center(polsys, sourceMask, Hrow);
MqGhost0 = local_charge_dipole_raw(polsys, sourceMask);

fprintf('\nGhost correction center:\n');
fprintf('  center = [%+.8e %+.8e %+.8e] bohr\n', ghostCenter);
fprintf('  raw Mq = [%+.8e %+.8e %+.8e] au | |Mq| = %.8e\n', ...
    MqGhost0, norm(MqGhost0));

%% ------------------------------------------------------------------------
% Mode comparison

fprintf('\n============================================================\n');
fprintf('Eext MODE / BOUNDARY COMPARISON\n');
fprintf('============================================================\n');

rows = {};
caseOutputs = struct([]);

for icase = 1:numel(cases)
    c = cases(icase);

    for ia = 1:numel(cfg.alphaList)
        alpha = cfg.alphaList(ia);

        fprintf('\n============================================================\n');
        fprintf('CASE %d/%d: %s\n', icase, numel(cases), c.label);
        fprintf('============================================================\n');
        fprintf('  Eext mode = %s\n', c.eextMode);
        fprintf('  boundary  = %s\n', c.boundary);
        fprintf('  alpha     = %.6f\n', alpha);
        fprintf('  rcut      = %.6f bohr\n', cfg.rcut);
        fprintf('  kcut      = %.6f bohr^-1\n', cfg.kcut);

        ewaldParams = struct();
        ewaldParams.alpha = alpha;
        ewaldParams.rcut = cfg.rcut;
        ewaldParams.kcut = cfg.kcut;
        ewaldParams.boundary = c.boundary;

        tCase = tic;

        %% ----------------------------------------------------------------
        % External field

        fprintf('[Eext] building %s external field...\n', c.eextMode);
        tE = tic;

        switch lower(c.eextMode)
            case 'nonperiodic'
                fieldNP = struct();
                fieldNP.exclude_self = true;
                fieldNP.use_thole_damping = cfg.use_thole;
                fieldNP.target_mask = targetMask;
                fieldNP.source_mask = sourceMask;

                Eext = thole.induced_field_from_charges(polsys, fieldNP);
                parts = local_make_nonperiodic_parts(polsys, Eext, targetMask, sourceMask);

            case 'periodic_ewald'
                fieldPer = struct();
                fieldPer.mode = 'periodic';
                fieldPer.exclude_self = true;
                fieldPer.use_thole_damping = cfg.use_thole;
                fieldPer.target_mask = targetMask;
                fieldPer.source_mask = sourceMask;
                fieldPer.verbose = false;
                fieldPer.kspace_mode = cfg.field_kspace_mode;
                fieldPer.k_block_size = cfg.field_k_block_size;
                fieldPer.kspace_memory_limit_gb = cfg.field_kspace_memory_limit_gb;
                fieldPer.ewald = ewaldParams;

                [Eext, parts] = thole.induced_field_from_charges_periodic(polsys, fieldPer);

            case 'periodic_p3m'
                p3mOpts = struct();
                p3mOpts.ewald = ewaldParams;
                p3mOpts.mesh_size = cfg.p3m.mesh_size;
                p3mOpts.assignment_order = cfg.p3m.assignment_order;
                p3mOpts.target_mask = targetMask;
                p3mOpts.source_mask = sourceMask;
                p3mOpts.exclude_self = true;
                p3mOpts.use_thole_damping = cfg.use_thole;
                p3mOpts.realspace_backend = cfg.p3m.realspace_backend;
                p3mOpts.deconvolve_assignment = cfg.p3m.deconvolve_assignment;
                p3mOpts.deconvolution_floor = cfg.p3m.deconvolution_floor;
                p3mOpts.derivative_mode = 'finite_difference';
                p3mOpts.fd_stencil = 'central2';
                p3mOpts.influence_mode = 'optimized';
                p3mOpts.alias_range = 2;
                p3mOpts.k_block_size = cfg.field_k_block_size;
                p3mOpts.kspace_memory_limit_gb = cfg.field_kspace_memory_limit_gb;
                p3mOpts.verbose = cfg.p3m.verbose;

                [Eext, parts] = p3m.compute_external_field_charges(polsys, p3mOpts);

            otherwise
                error('Unknown eextMode: %s', c.eextMode);
        end

        timeEext = toc(tE);

        fprintf('  Eext done in %.6f s\n', timeEext);
        fprintf('  ||Eext||_F    = %.16e\n', norm(Eext, 'fro'));
        fprintf('  ||Ereal||_F   = %.16e\n', norm(parts.real, 'fro'));
        fprintf('  ||Erecip||_F  = %.16e\n', norm(parts.recip, 'fro'));
        fprintf('  ||Esurf_q||_F = %.16e\n', norm(parts.surf, 'fro'));

        if isfield(parts, 'Mq')
            fprintf('  Mq(parts) = [%+.8e %+.8e %+.8e] au\n', parts.Mq);
        end
        if isfield(parts, 'Esurf_q')
            fprintf('  Esurf_q   = [%+.8e %+.8e %+.8e] au | norm = %.16e au\n', ...
                parts.Esurf_q, norm(parts.Esurf_q));
        end

        %% ----------------------------------------------------------------
        % Prepare problem

        sorParams = cfg.sor;
        sorParams.use_thole = cfg.use_thole;
        sorParams.kspace_mode = cfg.sor.kspace_mode;
        sorParams.k_block_size = cfg.sor.k_block_size;
        sorParams.kspace_memory_limit_gb = cfg.sor.kspace_memory_limit_gb;
        sorParams.use_mex_kspace = cfg.sor.use_mex_kspace;

        fprintf('[SOR] preparing SCF problem...\n');
        tPrep = tic;
        problem = thole.prepare_scf_problem(polsys, Eext, sorParams);
        timeProblemPrep = toc(tPrep);

        fprintf('  problem prep done in %.6f s | nPolSites = %d\n', ...
            timeProblemPrep, problem.nPolSites);

        %% ----------------------------------------------------------------
        % Prebuild periodic real-space dipole row cache with MEX

        fprintf('[SOR] building periodic real-space dipole row cache...\n');
        tRow = tic;

        rowOpts = struct();
        rowOpts.profile = true;
        rowOpts.use_mex = cfg.sor.use_mex_realspace;
        rowOpts.use_thole = cfg.use_thole;

        rowCache = geom.build_active_row_cache_periodic( ...
            polsys, problem, ewaldParams, rowOpts);

        timeRowCache = toc(tRow);

        nRealDipole = local_row_cache_count(rowCache);

        fprintf('  real-space row cache done in %.6f s\n', timeRowCache);
        fprintf('  nReal_dipole = %d\n', nRealDipole);

        sorParams.realspace_dipole_row_cache = rowCache;

        %% ----------------------------------------------------------------
        % Periodic matrix-free SOR

        fprintf('[SOR] solving matrix-free periodic SOR...\n');
        tSOR = tic;
        [mu, scf] = thole.solve_scf_iterative_periodic_sor( ...
            polsys, Eext, ewaldParams, sorParams);
        timeSOR = toc(tSOR);

        %% ----------------------------------------------------------------
        % Diagnostics

        polMask = logical(polsys.site_is_polarizable(:));
        muPol = mu(polMask, :);
        EextPol = Eext(polMask, :);

        E_stationary_Ha = -0.5 * sum(sum(muPol .* EextPol));
        E_stationary_eV = E_stationary_Ha * 27.211386245988;

        muMag = vecnorm(muPol, 2, 2);
        EextMag = vecnorm(EextPol, 2, 2);
        siteEnergyEV = -0.5 * sum(muPol .* EextPol, 2) * 27.211386245988;

        rmsMuAU = sqrt(mean(muMag.^2));
        rmsMuDebye = rmsMuAU * 2.541746473;

        maxMuAU = max(muMag);
        maxMuDebye = maxMuAU * 2.541746473;

        rmsEext = sqrt(mean(EextMag.^2));
        maxEext = max(EextMag);

        meanSiteEnergyEV = mean(siteEnergyEV);
        minSiteEnergyEV = min(siteEnergyEV);
        maxSiteEnergyEV = max(siteEnergyEV);

        [muSorted, muOrder] = sort(muMag, 'descend');
        nTopMu = min(cfg.nTopSites, numel(muSorted));
        topLocal = muOrder(1:nTopMu);
        polGlobalSites = find(polMask);
        topGlobal = polGlobalSites(topLocal);

        [Tmol, nTopMol] = local_molecule_energy_table(polsys, polMask, siteEnergyEV, cfg.nTopMolecules);

        surf = local_surface_energy_diagnostic(polsys, parts, muPol, c.boundary);

        %% ----------------------------------------------------------------
        % Point-dipole ghost-image corrections

        MqGhost = local_charge_dipole_raw(polsys, sourceMask);
        MmuGhost = sum(muPol, 1);
        MtotGhost = MqGhost + MmuGhost;

        ghostParams = struct();

        if cfg.ghost.use_case_ewald_params
            ghostParams.alpha = ewaldParams.alpha;
            ghostParams.rcut = ewaldParams.rcut;
            ghostParams.kcut = ewaldParams.kcut;
            ghostParams.boundary = ewaldParams.boundary;
        else
            ghostParams.alpha = cfg.ghost.alpha;
            ghostParams.rcut = cfg.ghost.rcut;
            ghostParams.kcut = cfg.ghost.kcut;
            ghostParams.boundary = ewaldParams.boundary;
        end

        ghostMu = local_point_dipole_ghost_energy_correction( ...
            polsys, polMask, muPol, ghostCenter, MmuGhost, Hrow, ghostParams);

        ghostTot = local_point_dipole_ghost_energy_correction( ...
            polsys, polMask, muPol, ghostCenter, MtotGhost, Hrow, ghostParams);

        E_ghostMu_eV  = E_stationary_eV + ghostMu.energyCorrectionEV;
        E_ghostTot_eV = E_stationary_eV + ghostTot.energyCorrectionEV;

        %% ----------------------------------------------------------------
        % Print result summary

        fprintf('\nResult summary for case: %s\n', c.label);
        fprintf('  boundary actually passed to SOR = %s\n', ewaldParams.boundary);
        fprintf('  converged = %d | nIter = %d | relres = %.16e\n', ...
            scf.converged, scf.nIter, scf.relres);
        fprintf('  ||Eext||_F = %.16e\n', norm(Eext, 'fro'));
        fprintf('  ||mu||_F   = %.16e\n', norm(mu, 'fro'));
        fprintf('  RMS |mu_i| = %.16e au = %.8f D\n', rmsMuAU, rmsMuDebye);
        fprintf('  max |mu_i| = %.16e au = %.8f D\n', maxMuAU, maxMuDebye);
        fprintf('  RMS |Eext_i| = %.16e au\n', rmsEext);
        fprintf('  max |Eext_i| = %.16e au\n', maxEext);
        fprintf('  E_stationary = %+0.16e Ha = %+0.8f eV\n', ...
            E_stationary_Ha, E_stationary_eV);
        fprintf('  mean site energy = %+0.8e eV/site\n', meanSiteEnergyEV);
        fprintf('  min/max site energy = %+0.8e / %+0.8e eV\n', ...
            minSiteEnergyEV, maxSiteEnergyEV);

        fprintf('  surface dipoles:\n');
        fprintf('    Mq(parts) = [%+.8e %+.8e %+.8e]\n', surf.Mq);
        fprintf('    Mmu       = [%+.8e %+.8e %+.8e]\n', surf.Mmu);
        fprintf('    Mtot      = [%+.8e %+.8e %+.8e]\n', surf.Mtot);
        fprintf('  surface energy diagnostic for boundary=%s:\n', c.boundary);
        fprintf('    Uqq    = %+0.8f eV\n', surf.Uqq_eV);
        fprintf('    Ucross = %+0.8f eV\n', surf.Ucross_eV);
        fprintf('    Umumu  = %+0.8f eV\n', surf.Umumu_eV);
        fprintf('    Utotal = %+0.8f eV\n', surf.Utotal_eV);

        fprintf('  point-dipole ghost-image corrections:\n');
        fprintf('    ghost center COM = [%+.8e %+.8e %+.8e] bohr\n', ghostCenter);
        fprintf('    Mq(raw)          = [%+.8e %+.8e %+.8e] au | |Mq| = %.8e\n', ...
            MqGhost, norm(MqGhost));
        fprintf('    Mmu              = [%+.8e %+.8e %+.8e] au | |Mmu| = %.8e\n', ...
            MmuGhost, norm(MmuGhost));
        fprintf('    Mtot             = [%+.8e %+.8e %+.8e] au | |Mtot| = %.8e\n', ...
            MtotGhost, norm(MtotGhost));

        fprintf('    induced-only image correction:\n');
        fprintf('      ||Eimg(Mmu)||_F = %.16e\n', norm(ghostMu.Eimg, 'fro'));
        fprintf('      max |Eimg_i|    = %.16e\n', max(vecnorm(ghostMu.Eimg(polMask,:), 2, 2)));
        fprintf('      dE_mu_img       = %+0.8f eV\n', ghostMu.energyCorrectionEV);
        fprintf('      E + dE_mu_img   = %+0.8f eV\n', E_ghostMu_eV);

        fprintf('    total-cell-dipole image correction:\n');
        fprintf('      ||Eimg(Mtot)||_F = %.16e\n', norm(ghostTot.Eimg, 'fro'));
        fprintf('      max |Eimg_i|     = %.16e\n', max(vecnorm(ghostTot.Eimg(polMask,:), 2, 2)));
        fprintf('      dE_tot_img       = %+0.8f eV\n', ghostTot.energyCorrectionEV);
        fprintf('      E + dE_tot_img   = %+0.8f eV\n', E_ghostTot_eV);

        if cfg.printTopSites
            fprintf('  top induced dipoles:\n');
            for kk = 1:nTopMu
                fprintf('    site %6d | |mu| = %.8e au = %.6f D | |Eext| = %.8e au | siteE = %+0.8e eV\n', ...
                    topGlobal(kk), muSorted(kk), muSorted(kk) * 2.541746473, ...
                    EextMag(topLocal(kk)), siteEnergyEV(topLocal(kk)));
            end
        end

        if cfg.printTopMolecules
            if nTopMol > 0
                fprintf('  top molecule energy contributions:\n');
                for kk = 1:nTopMol
                    fprintf('    mol %6d | E_mol = %+0.8f eV | abs = %.8f eV\n', ...
                        Tmol.molID(kk), Tmol.energyEV(kk), Tmol.absEnergyEV(kk));
                end
                fprintf('  sum molecule energies = %+0.8f eV\n', sum(Tmol.energyEV));
            else
                fprintf('  molecule energy diagnostics unavailable: polsys.site_mol_id missing.\n');
            end
        end

        timeTotal = toc(tCase);

        rows(end+1,:) = { ...
            c.label, c.eextMode, c.boundary, alpha, cfg.rcut, cfg.kcut, ...
            norm(Eext,'fro'), norm(parts.real,'fro'), norm(parts.recip,'fro'), norm(parts.surf,'fro'), ...
            local_get_field_or_nan(parts, 'nRealEntries'), local_get_field_or_nan(parts, 'nK'), ...
            nRealDipole, ...
            scf.converged, scf.nIter, scf.relres, ...
            norm(mu,'fro'), rmsMuAU, rmsMuDebye, maxMuAU, maxMuDebye, ...
            rmsEext, maxEext, ...
            E_stationary_Ha, E_stationary_eV, ...
            meanSiteEnergyEV, minSiteEnergyEV, maxSiteEnergyEV, ...
            surf.Uqq_eV, surf.Ucross_eV, surf.Umumu_eV, surf.Utotal_eV, ...
            norm(MqGhost), norm(MmuGhost), norm(MtotGhost), ...
            ghostMu.energyCorrectionEV, E_ghostMu_eV, ...
            ghostTot.energyCorrectionEV, E_ghostTot_eV, ...
            norm(ghostMu.Eimg,'fro'), max(vecnorm(ghostMu.Eimg(polMask,:),2,2)), ...
            norm(ghostTot.Eimg,'fro'), max(vecnorm(ghostTot.Eimg(polMask,:),2,2)), ...
            local_get_part_time(parts, 'time_real'), local_get_part_time(parts, 'time_recip'), local_get_part_time(parts, 'time_mesh'), ...
            timeEext, timeProblemPrep, timeRowCache, timeSOR, timeTotal}; %#ok<SAGROW>

        out = struct();
        out.case = c;
        out.alpha = alpha;
        out.Eext = Eext;
        out.parts = parts;
        out.problem = problem;
        out.rowCache = rowCache;
        out.mu = mu;
        out.scf = scf;
        out.E_stationary_eV = E_stationary_eV;
        out.surface = surf;
        out.ghostCenter = ghostCenter;
        out.MqGhost = MqGhost;
        out.MmuGhost = MmuGhost;
        out.MtotGhost = MtotGhost;
        out.ghostMu = ghostMu;
        out.ghostTot = ghostTot;
        out.E_ghostMu_eV = E_ghostMu_eV;
        out.E_ghostTot_eV = E_ghostTot_eV;

        caseOutputs(end+1).out = out; %#ok<SAGROW>
    end
end

summary = cell2table(rows, 'VariableNames', { ...
    'caseLabel', 'eextMode', 'boundary', 'alpha', 'rcut', 'kcut', ...
    'normEext', 'normEreal', 'normErecip', 'normEsurfQ', ...
    'nRealField', 'nKField', ...
    'nRealDipole', ...
    'converged', 'nIter', 'relres', ...
    'normMu', 'rmsMuAU', 'rmsMuDebye', 'maxMuAU', 'maxMuDebye', ...
    'rmsEext', 'maxEext', ...
    'energyHa', 'energyEV', ...
    'meanSiteEnergyEV', 'minSiteEnergyEV', 'maxSiteEnergyEV', ...
    'UqqEV', 'UcrossEV', 'UmumuEV', 'UtotalEV', ...
    'normMqGhost', 'normMmuGhost', 'normMtotGhost', ...
    'dEghostMuEV', 'energyGhostMuEV', ...
    'dEghostTotEV', 'energyGhostTotEV', ...
    'normEimgMu', 'maxEimgMu', ...
    'normEimgTot', 'maxEimgTot', ...
    'timeEextReal', 'timeEextRecip', 'timeEextMesh', ...
    'timeEext', 'timeProblemPrep', 'timeRowCache', 'timeSOR', 'timeTotal'});

fprintf('\n============================================================\n');
fprintf('MODE COMPARISON SUMMARY\n');
fprintf('============================================================\n');
disp(summary);

fprintf('\nCompact comparison:\n');
for i = 1:height(summary)
    fprintf(['  %-45s | Eext %-13s | boundary %-8s | alpha %.3f | ' ...
             '||Eext|| %.8e | E %.8f eV | E+dE_mu %.8f eV | E+dE_tot %.8f eV | ' ...
             '||mu|| %.8e | RMSmu %.6f D | relres %.3e | iter %d | Eext time %.3f s\n'], ...
        summary.caseLabel{i}, summary.eextMode{i}, summary.boundary{i}, summary.alpha(i), ...
        summary.normEext(i), summary.energyEV(i), summary.energyGhostMuEV(i), summary.energyGhostTotEV(i), ...
        summary.normMu(i), summary.rmsMuDebye(i), ...
        summary.relres(i), summary.nIter(i), summary.timeEext(i));
end

fprintf('\nExternal-field comparison notes:\n');
idxEwaldTin = find(strcmp(summary.eextMode, 'periodic_ewald') & strcmp(summary.boundary, 'tinfoil'), 1);
idxP3MTin = find(strcmp(summary.eextMode, 'periodic_p3m') & strcmp(summary.boundary, 'tinfoil'), 1);
if ~isempty(idxEwaldTin) && ~isempty(idxP3MTin)
    fprintf('  periodic Ewald/tinfoil Eext time = %.6f s\n', summary.timeEext(idxEwaldTin));
    fprintf('  periodic P3M/tinfoil   Eext time = %.6f s\n', summary.timeEext(idxP3MTin));
    fprintf('  speedup Ewald/P3M = %.6f x\n', summary.timeEext(idxEwaldTin) / summary.timeEext(idxP3MTin));

    fprintf('  |E_P3M - E_Ewald| energy difference = %.8e eV\n', ...
        abs(summary.energyEV(idxP3MTin) - summary.energyEV(idxEwaldTin)));
    fprintf('  |E_P3M+dE_mu - E_Ewald+dE_mu| = %.8e eV\n', ...
        abs(summary.energyGhostMuEV(idxP3MTin) - summary.energyGhostMuEV(idxEwaldTin)));
end

fprintf('\nEnergy ratios relative to first raw case:\n');
E0 = summary.energyEV(1);
for i = 1:height(summary)
    fprintf('  %-45s | E/E_case1 = %.8f | E_muCorr/E_case1 = %.8f | E_totCorr/E_case1 = %.8f\n', ...
        summary.caseLabel{i}, ...
        summary.energyEV(i) / E0, ...
        summary.energyGhostMuEV(i) / E0, ...
        summary.energyGhostTotEV(i) / E0);
end

fprintf('\nRecommended interpretation:\n');
fprintf('  - For nonperiodic_Eext__periodic_SOR_tinfoil, use E+dE_mu as the ghost-corrected hybrid estimate.\n');
fprintf('  - The periodic_P3M_Eext row should match periodic_Ewald_Eext/tinfoil, but with faster Eext construction.\n');
fprintf('  - Periodic-source rows remain periodic CT-array diagnostics, not isolated-CT production estimates.\n');

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

function parts = local_make_nonperiodic_parts(sys, Eext, targetMask, sourceMask)
    nSites = size(sys.site_pos, 1);

    parts = struct();
    parts.real = Eext;
    parts.recip = zeros(nSites, 3);
    parts.surf = zeros(nSites, 3);
    parts.total = Eext;

    parts.Mq = [0.0, 0.0, 0.0];
    parts.Esurf_q = [0.0, 0.0, 0.0];
    parts.surf_coeff = 0.0;

    parts.target_mask = targetMask;
    parts.source_mask = sourceMask;
    parts.target_sites = find(targetMask);
    parts.source_sites = find(sourceMask);
    parts.qtot = sum(sys.site_charge(sourceMask));

    parts.nTargets = nnz(targetMask);
    parts.nSources = nnz(sourceMask);
    parts.nRealEntries = NaN;
    parts.nK = 0;

    parts.alpha = NaN;
    parts.rcut = NaN;
    parts.kcut = NaN;
    parts.boundary = 'nonperiodic';
    parts.use_thole_damping = NaN;
    parts.exclude_self = true;
    parts.real_only = false;

    parts.storage_mode = 'none';
    parts.kspace_mode_requested = 'none';
    parts.k_block_size = NaN;
    parts.kspace_memory_limit_gb = NaN;
    parts.estimated_full_gb = 0.0;

    parts.time_real = NaN;
    parts.time_recip = 0.0;
    parts.time_mesh = NaN;
    parts.time_total = NaN;

    parts.realCache = struct();
end

function n = local_row_cache_count(rowCache)
    if isfield(rowCache, 'nInteractions') && ~isempty(rowCache.nInteractions)
        n = rowCache.nInteractions;
    elseif isfield(rowCache, 'nEntriesDirected') && ~isempty(rowCache.nEntriesDirected)
        n = rowCache.nEntriesDirected;
    elseif isfield(rowCache, 'nEntries') && ~isempty(rowCache.nEntries)
        n = rowCache.nEntries;
    else
        n = NaN;
    end
end

function [Tmol, nTopMol] = local_molecule_energy_table(polsys, polMask, siteEnergyEV, nTop)
    if isfield(polsys, 'site_mol_id') && ~isempty(polsys.site_mol_id)
        molIDPol = polsys.site_mol_id(polMask);

        [Gmol, molIDs] = findgroups(molIDPol);
        molEnergyEV = splitapply(@sum, siteEnergyEV, Gmol);
        molAbsEnergyEV = abs(molEnergyEV);

        Tmol = table(molIDs, molEnergyEV, molAbsEnergyEV, ...
            'VariableNames', {'molID', 'energyEV', 'absEnergyEV'});

        Tmol = sortrows(Tmol, 'absEnergyEV', 'descend');
        nTopMol = min(nTop, height(Tmol));
    else
        Tmol = table();
        nTopMol = 0;
    end
end

function surf = local_surface_energy_diagnostic(polsys, parts, muPol, boundary)
    Hrow = local_get_direct_lattice(polsys);
    Hcol = Hrow.';
    V = abs(det(Hcol));

    Mq = [0.0, 0.0, 0.0];
    if isfield(parts, 'Mq') && ~isempty(parts.Mq)
        Mq = parts.Mq;
    end

    Mmu = sum(muPol, 1);
    Mtot = Mq + Mmu;

    switch lower(boundary)
        case 'tinfoil'
            surf_energy_coeff = 0.0;

        case 'vacuum'
            surf_energy_coeff = 2*pi/(3*V);

        otherwise
            error('Unknown boundary: %s', boundary);
    end

    Uqq_Ha = surf_energy_coeff * dot(Mq, Mq);
    Umumu_Ha = surf_energy_coeff * dot(Mmu, Mmu);
    Ucross_Ha = 2 * surf_energy_coeff * dot(Mq, Mmu);
    Utotal_Ha = surf_energy_coeff * dot(Mtot, Mtot);

    Ha_to_eV = 27.211386245988;

    surf = struct();
    surf.Mq = Mq;
    surf.Mmu = Mmu;
    surf.Mtot = Mtot;
    surf.Uqq_Ha = Uqq_Ha;
    surf.Umumu_Ha = Umumu_Ha;
    surf.Ucross_Ha = Ucross_Ha;
    surf.Utotal_Ha = Utotal_Ha;
    surf.Uqq_eV = Uqq_Ha * Ha_to_eV;
    surf.Umumu_eV = Umumu_Ha * Ha_to_eV;
    surf.Ucross_eV = Ucross_Ha * Ha_to_eV;
    surf.Utotal_eV = Utotal_Ha * Ha_to_eV;
end

function center = local_charged_pair_com_center(sys, sourceMask, Hrow)
%LOCAL_CHARGED_PAIR_COM_CENTER Coherent COM midpoint of charged +/- sites.
%
% Uses source charge sign to split the two charged molecules/sites. The
% negative-site COM is minimum-imaged relative to the positive-site COM, then
% the center is the midpoint.
%
% This center is only the position at which we place the ghost point dipole.
% The CT dipole Mq used in diagnostics is computed by local_charge_dipole_raw.

    pos = sys.site_pos;
    q = sys.site_charge(:);

    src = find(sourceMask);
    qsrc = q(src);
    posSrc = pos(src, :);

    plus = qsrc > 0;
    minus = qsrc < 0;

    if ~any(plus) || ~any(minus)
        error('local_charged_pair_com_center:NeedPlusMinus', ...
            'Expected both positive and negative source charges.');
    end

    rPlus = mean(posSrc(plus, :), 1);
    rMinusRaw = mean(posSrc(minus, :), 1);

    dMinus = local_minimum_image_displacements(rMinusRaw, rPlus, Hrow);
    rMinus = rPlus + dMinus;

    center = 0.5 * (rPlus + rMinus);
end

function Mq = local_charge_dipole_raw(sys, sourceMask)
%LOCAL_CHARGE_DIPOLE_RAW Raw charge dipole from current Cartesian branch.
%
% For the chosen charged pair in this workflow, the raw Cartesian branch is
% expected to represent the physical same-stack CT separation.

    src = find(sourceMask);
    qsrc = sys.site_charge(src);
    rsrc = sys.site_pos(src, :);

    Mq = sum(qsrc .* rsrc, 1);
end

function ghost = local_point_dipole_ghost_energy_correction( ...
    sys, polMask, muPol, center, Mghost, Hrow, params)
%LOCAL_POINT_DIPOLE_GHOST_ENERGY_CORRECTION
%
% Computes image-only point-dipole field:
%
%   Eimg = E_periodic_point_dipole(Mghost) - E_central_bare(Mghost)
%
% using minimum-image displacements relative to center.
%
% Then returns:
%
%   dE = +1/2 sum_i mu_i . Eimg_i
%
% This is a post-hoc correction estimate for ghost-image contributions.

    targetSites = find(polMask);
    targetPos = sys.site_pos(targetSites, :);

    alpha = params.alpha;
    rcut = params.rcut;
    kcut = params.kcut;
    boundary = params.boundary;

    Hcol = Hrow.';

    dMin = local_minimum_image_displacements(targetPos, center, Hrow);

    targetLocal = dMin;
    centerLocal = [0.0, 0.0, 0.0];

    [EperTarget, parts] = local_periodic_point_dipole_field( ...
        targetLocal, centerLocal, Mghost, Hcol, alpha, rcut, kcut, boundary);

    EbareTarget = local_bare_point_dipole_field(targetLocal, centerLocal, Mghost);

    EimgTarget = EperTarget - EbareTarget;

    Eper = zeros(size(sys.site_pos));
    Ebare = zeros(size(sys.site_pos));
    Eimg = zeros(size(sys.site_pos));

    Eper(targetSites, :) = EperTarget;
    Ebare(targetSites, :) = EbareTarget;
    Eimg(targetSites, :) = EimgTarget;

    dE_Ha = +0.5 * sum(sum(muPol .* EimgTarget));
    dE_eV = dE_Ha * 27.211386245988;

    ghost = struct();
    ghost.center = center;
    ghost.Mghost = Mghost;
    ghost.alpha = alpha;
    ghost.rcut = rcut;
    ghost.kcut = kcut;
    ghost.boundary = boundary;
    ghost.Eper = Eper;
    ghost.Ebare = Ebare;
    ghost.Eimg = Eimg;
    ghost.parts = parts;
    ghost.energyCorrectionHa = dE_Ha;
    ghost.energyCorrectionEV = dE_eV;
end

function dMin = local_minimum_image_displacements(targetPos, center, Hrow)
%LOCAL_MINIMUM_IMAGE_DISPLACEMENTS Minimum-image displacement to center.
%
% Hrow stores direct lattice vectors as rows:
%
%   cart = frac * Hrow
%
% For each target, this returns the Cartesian displacement from center to
% target after wrapping the fractional displacement into [-1/2, 1/2).

    d = targetPos - center;

    frac = d / Hrow;
    frac = frac - round(frac);

    dMin = frac * Hrow;
end

function [Eper, parts] = local_periodic_point_dipole_field( ...
    targetPos, center, M, Hcol, alpha, rcut, kcut, boundary)

    [Ereal, realParts] = local_realspace_point_dipole_field( ...
        targetPos, center, M, Hcol, alpha, rcut);

    [Erecip, recipParts] = local_recip_point_dipole_field( ...
        targetPos, center, M, Hcol, alpha, kcut);

    Esurf = zeros(size(targetPos));

    V = abs(det(Hcol));
    switch lower(boundary)
        case 'tinfoil'
            % no-op

        case 'vacuum'
            Es = -(4*pi/(3*V)) * M;
            Esurf = repmat(Es, size(targetPos,1), 1);

        otherwise
            error('local_periodic_point_dipole_field:BadBoundary', ...
                'boundary must be ''tinfoil'' or ''vacuum''.');
    end

    Eper = Ereal + Erecip + Esurf;

    parts = struct();
    parts.real = Ereal;
    parts.recip = Erecip;
    parts.surf = Esurf;
    parts.realParts = realParts;
    parts.recipParts = recipParts;
end

function [Ereal, parts] = local_realspace_point_dipole_field( ...
    targetPos, center, M, Hcol, alpha, rcut)

    nT = size(targetPos, 1);
    Ereal = zeros(nT, 3);

    % Conservative bounds. This may include extra translations, which are
    % filtered by rcut.
    vecLens = vecnorm(Hcol, 2, 1);
    nMax = ceil(rcut / min(vecLens)) + 3;
    nRange = [nMax nMax nMax];

    nTerms = 0;

    for n1 = -nRange(1):nRange(1)
        for n2 = -nRange(2):nRange(2)
            for n3 = -nRange(3):nRange(3)
                R = (Hcol * [n1; n2; n3]).';

                dr = targetPos - (center + R);
                r2 = sum(dr.^2, 2);

                keep = (r2 > 1e-24) & (r2 <= rcut^2);
                if ~any(keep)
                    continue;
                end

                drk = dr(keep, :);
                r2k = r2(keep);
                rk = sqrt(r2k);

                ar = alpha .* rk;
                expTerm = exp(-(alpha^2) .* r2k);

                B = erfc(ar) ./ (rk.^3) + ...
                    (2 * alpha / sqrt(pi)) .* expTerm ./ r2k;

                C = 3 * erfc(ar) ./ (rk.^5) + ...
                    (2 * alpha / sqrt(pi)) .* expTerm .* ...
                    (3 ./ (rk.^4) + 2 * alpha^2 ./ r2k);

                coeff_iso = -B;
                coeff_dyad = C;

                muDotR = drk * M.';
                contrib = coeff_iso .* M + coeff_dyad .* (muDotR .* drk);

                Ereal(keep, :) = Ereal(keep, :) + contrib;

                nTerms = nTerms + nnz(keep);
            end
        end
    end

    parts = struct();
    parts.nTerms = nTerms;
end

function [Erecip, parts] = local_recip_point_dipole_field( ...
    targetPos, center, M, Hcol, alpha, kcut)

    nT = size(targetPos, 1);
    Erecip = zeros(nT, 3);

    V = abs(det(Hcol));

    [kvecs, meta] = ewald.enumerate_kvecs_triclinic(Hcol, kcut);
    nk = size(kvecs, 1);

    if nk == 0
        parts = struct('nK', 0);
        return;
    end

    k2 = meta.k2(:);

    % Half-k convention matching the periodic dipole operator:
    %
    %   T_recip(d) =
    %       -8*pi/V sum_half exp(-k^2/4a^2)/k^2
    %       cos(k.d) k k^T
    two_pref = -(8*pi/V) * exp(-k2 ./ (4 * alpha^2)) ./ k2;

    d = targetPos - center;
    phase = d * kvecs.';
    cosPhase = cos(phase);

    Mk = M * kvecs.';              % 1 x nk
    W = cosPhase .* (two_pref.' .* Mk);

    Erecip(:,1) = W * kvecs(:,1);
    Erecip(:,2) = W * kvecs(:,2);
    Erecip(:,3) = W * kvecs(:,3);

    parts = struct();
    parts.nK = nk;
end

function Ebare = local_bare_point_dipole_field(targetPos, center, M)
    nT = size(targetPos, 1);
    Ebare = zeros(nT, 3);

    I3 = eye(3);

    dr = targetPos - center;
    r2 = sum(dr.^2, 2);

    for i = 1:nT
        if r2(i) <= 1e-24
            continue;
        end

        r = sqrt(r2(i));
        x = dr(i, :).';

        T = -(1 / r^3) * I3 + 3 * (x * x.') / r^5;
        Ebare(i, :) = (T * M.').';
    end
end

function val = local_get_field_or_nan(s, name)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = NaN;
    end
end

function val = local_get_part_time(parts, name)
    if isfield(parts, name) && ~isempty(parts.(name))
        val = parts.(name);
    elseif strcmp(name, 'time_mesh') && isfield(parts, 'time_recip')
        % Direct Ewald has time_recip; P3M has time_mesh.
        val = NaN;
    else
        val = NaN;
    end
end