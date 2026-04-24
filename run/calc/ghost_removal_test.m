%% test_point_dipole_image_ewald_real_geometry
% Real-geometry debug for induced-dipole ghost-image correction.
%
% Builds the real Polarize system, solves one SCF case to obtain mu, then
% approximates the periodically repeated induced cloud by one point dipole
%
%   Mmu = sum_i mu_i
%
% placed at a chosen center. It computes
%
%   E_img(r_i) = E_periodic_point_dipole(r_i; Mmu)
%              - E_central_bare_point_dipole(r_i; Mmu)
%
% on the real polarizable sites.
%
% Critical correction in this version:
%   The point-dipole field is evaluated using minimum-image displacements
%   from the chosen center to each target site. This avoids mixing the
%   periodic point-dipole Ewald field with arbitrary unwrapped/display
%   coordinates.
%
% Diagnostic only. Do not use the correction until:
%   - E_img is stable with ghost alpha / rcut / kcut
%   - E_img has a plausible magnitude
%   - the center convention is sane

clear; clc; close all;

fprintf('=== real-geometry point-dipole image Ewald debug ===\n');

%% ------------------------------------------------------------------------
% User controls

cfg = struct();

cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename   = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

cfg.supercellSize = [2 5 1];
cfg.bondScale     = 1.20;
cfg.pairCharges   = [+1 -1];

cfg.use_thole = true;
cfg.verbose = true;

% SCF case used to generate the "relevant" induced dipole.
% Use the hybrid case first because this is closest to the isolated-source
% diagnostic we have been discussing.
cfg.scfEextMode = 'nonperiodic';   % 'nonperiodic' or 'periodic'
cfg.scfBoundary = 'tinfoil';       % 'tinfoil' or 'vacuum'

cfg.scfAlpha = 0.35;
cfg.scfRcut  = 11.5;
cfg.scfKcut  = 3.5;

cfg.field_kspace_mode = 'full';
cfg.field_k_block_size = 2048;
cfg.field_kspace_memory_limit_gb = 8;

cfg.sor = struct();
cfg.sor.use_thole = cfg.use_thole;
cfg.sor.verbose = true;
cfg.sor.printEvery = 10;
cfg.sor.residualEvery = 10;
cfg.sor.tol = 1e-6;
cfg.sor.maxIter = 500;
cfg.sor.omega = 0.55;
cfg.sor.stopMetric = 'max_dmu';

cfg.sor.use_mex_realspace = true;

cfg.sor.kspace_mode = 'full';
cfg.sor.k_block_size = 2048;
cfg.sor.kspace_memory_limit_gb = 8;
cfg.sor.use_mex_kspace = false;

% Ghost-image debug sweeps.
%
% These are for the point-dipole ghost correction only.
% They do not need to equal the atomistic SOR row-cache rcut.
cfg.ghostAlphaList = [0.20 0.25 0.30 0.35 0.40 0.45];
cfg.ghostRcutList  = [11.5 20.0 40.0 60.0];
cfg.ghostKcutList  = [3.5 5.0 7.0];

% Main alpha sweep uses the largest chosen rcut/kcut by default.
cfg.ghostAlphaSweepRcut = cfg.ghostRcutList(end);
cfg.ghostAlphaSweepKcut = cfg.ghostKcutList(end);

% Center choices to test.
% charged_mean: mean position of all charged sites.
% cell_center : geometric center of the supercell parallelepiped.
cfg.centerModes = {'charged_mean', 'cell_center'};

% Print largest image fields.
cfg.nTopTargets = 12;

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
% Build system

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
% Choose charged pair

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
% Apply charges and extract polarization system

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

targetMask = logical(polsys.site_is_polarizable(:));
sourceMask = abs(polsys.site_charge(:)) > 0;

fprintf('\nCanonical polarization system:\n');
fprintf('  n_sites             = %d\n', polsys.n_sites);
fprintf('  n_polarizable_sites = %d\n', nnz(targetMask));
fprintf('  n_source_sites      = %d\n', nnz(sourceMask));
fprintf('  total source charge = %+0.16e\n', sum(polsys.site_charge(sourceMask)));

if nnz(targetMask & sourceMask) ~= 0
    error('Charged source sites are still polarizable.');
end
if abs(sum(polsys.site_charge(sourceMask))) > 1e-10
    error('Selected source charges are not neutral.');
end

Hrow = local_get_direct_lattice(polsys);
Hcol = Hrow.';
V = abs(det(Hcol));
Lmin = geom.shortest_lattice_translation(Hcol);

fprintf('\nLattice sanity:\n');
fprintf('  volume = %.16e bohr^3\n', V);
fprintf('  Lmin   = %.16e bohr\n', Lmin);
fprintf('  Lmin/2 = %.16e bohr\n', 0.5 * Lmin);

%% ------------------------------------------------------------------------
% Build Eext for the selected SCF case

fprintf('\n============================================================\n');
fprintf('BUILDING SCF CASE FOR RELEVANT INDUCED DIPOLES\n');
fprintf('============================================================\n');

ewaldParams = struct();
ewaldParams.alpha = cfg.scfAlpha;
ewaldParams.rcut = cfg.scfRcut;
ewaldParams.kcut = cfg.scfKcut;
ewaldParams.boundary = cfg.scfBoundary;

fprintf('  Eext mode = %s\n', cfg.scfEextMode);
fprintf('  boundary  = %s\n', cfg.scfBoundary);
fprintf('  alpha     = %.6f\n', cfg.scfAlpha);
fprintf('  rcut      = %.6f bohr\n', cfg.scfRcut);
fprintf('  kcut      = %.6f bohr^-1\n', cfg.scfKcut);

switch lower(cfg.scfEextMode)
    case 'nonperiodic'
        fieldNP = struct();
        fieldNP.exclude_self = true;
        fieldNP.use_thole_damping = cfg.use_thole;
        fieldNP.target_mask = targetMask;
        fieldNP.source_mask = sourceMask;

        Eext = thole.induced_field_from_charges(polsys, fieldNP);

        parts = local_make_nonperiodic_parts(polsys, Eext, targetMask, sourceMask);

    case 'periodic'
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

    otherwise
        error('Unknown cfg.scfEextMode: %s', cfg.scfEextMode);
end

fprintf('\nEext summary:\n');
fprintf('  ||Eext||_F    = %.16e\n', norm(Eext, 'fro'));
fprintf('  ||Ereal||_F   = %.16e\n', norm(parts.real, 'fro'));
fprintf('  ||Erecip||_F  = %.16e\n', norm(parts.recip, 'fro'));
fprintf('  ||Esurf_q||_F = %.16e\n', norm(parts.surf, 'fro'));

%% ------------------------------------------------------------------------
% Solve periodic SOR to get mu

fprintf('\n[SOR] preparing problem and row cache...\n');

sorParams = cfg.sor;
sorParams.use_thole = cfg.use_thole;
sorParams.kspace_mode = cfg.sor.kspace_mode;
sorParams.k_block_size = cfg.sor.k_block_size;
sorParams.kspace_memory_limit_gb = cfg.sor.kspace_memory_limit_gb;
sorParams.use_mex_kspace = cfg.sor.use_mex_kspace;

problem = thole.prepare_scf_problem(polsys, Eext, sorParams);

rowOpts = struct();
rowOpts.profile = true;
rowOpts.use_mex = cfg.sor.use_mex_realspace;
rowOpts.use_thole = cfg.use_thole;

rowCache = geom.build_active_row_cache_periodic( ...
    polsys, problem, ewaldParams, rowOpts);

sorParams.realspace_dipole_row_cache = rowCache;

fprintf('\n[SOR] solving...\n');
[mu, scf] = thole.solve_scf_iterative_periodic_sor( ...
    polsys, Eext, ewaldParams, sorParams);

polMask = targetMask;
muPol = mu(polMask, :);
EextPol = Eext(polMask, :);

E_stationary_Ha = -0.5 * sum(sum(muPol .* EextPol));
E_stationary_eV = E_stationary_Ha * 27.211386245988;

Mmu = sum(muPol, 1);
muMag = vecnorm(muPol, 2, 2);

fprintf('\nSCF result:\n');
fprintf('  converged = %d | nIter = %d | relres = %.16e\n', ...
    scf.converged, scf.nIter, scf.relres);
fprintf('  E_stationary = %+0.16e Ha = %+0.8f eV\n', ...
    E_stationary_Ha, E_stationary_eV);
fprintf('  ||mu||_F = %.16e\n', norm(mu, 'fro'));
fprintf('  Mmu = [%+.8e %+.8e %+.8e] au\n', Mmu);
fprintf('  |Mmu| = %.16e au\n', norm(Mmu));
fprintf('  RMS |mu_i| = %.16e au = %.8f D\n', ...
    sqrt(mean(muMag.^2)), sqrt(mean(muMag.^2))*2.541746473);

%% ------------------------------------------------------------------------
% Real-geometry point-dipole image tests

fprintf('\n============================================================\n');
fprintf('REAL-GEOMETRY POINT-DIPOLE IMAGE TESTS\n');
fprintf('============================================================\n');

targetSites = find(polMask);
targetPos = polsys.site_pos(targetSites, :);

for ic = 1:numel(cfg.centerModes)
    centerMode = cfg.centerModes{ic};
    center = local_choose_center(polsys, sourceMask, Hrow, centerMode);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('CENTER MODE: %s\n', centerMode);
    fprintf('------------------------------------------------------------\n');
    fprintf('  center = [%+.8e %+.8e %+.8e] bohr\n', center);

    dRaw = targetPos - center;
    dMin = local_minimum_image_displacements(targetPos, center, Hrow);

    dRawMag = vecnorm(dRaw, 2, 2);
    dMinMag = vecnorm(dMin, 2, 2);

    fprintf('  raw target distance to center:\n');
    fprintf('    min    = %.8e bohr\n', min(dRawMag));
    fprintf('    median = %.8e bohr\n', median(dRawMag));
    fprintf('    max    = %.8e bohr\n', max(dRawMag));

    fprintf('  minimum-image target distance to center:\n');
    fprintf('    min    = %.8e bohr\n', min(dMinMag));
    fprintf('    median = %.8e bohr\n', median(dMinMag));
    fprintf('    max    = %.8e bohr\n', max(dMinMag));

    [dSorted, dOrder] = sort(dMinMag, 'ascend');
    fprintf('  nearest minimum-image targets to center:\n');
    for kk = 1:min(cfg.nTopTargets, numel(dSorted))
        fprintf('    site %6d | dist_min_image = %.8e bohr | dist_raw = %.8e bohr\n', ...
            targetSites(dOrder(kk)), dSorted(kk), dRawMag(dOrder(kk)));
    end

    %% --------------------------------------------------------------------
    % Rcut/kcut convergence at fixed alpha

    alpha0 = cfg.scfAlpha;

    rows = {};

    fprintf('\n  Rcut/Kcut sweep at alpha = %.6f:\n', alpha0);

    for ir = 1:numel(cfg.ghostRcutList)
        rcutGhost = cfg.ghostRcutList(ir);

        for ik = 1:numel(cfg.ghostKcutList)
            kcutGhost = cfg.ghostKcutList(ik);

            ghost = local_point_dipole_image_field_real_geometry( ...
                targetPos, center, Mmu, Hrow, alpha0, rcutGhost, kcutGhost, cfg.scfBoundary);

            dE_Ha = +0.5 * sum(sum(muPol .* ghost.Eimg));
            dE_eV = dE_Ha * 27.211386245988;

            rows(end+1,:) = { ...
                centerMode, alpha0, rcutGhost, kcutGhost, ...
                norm(ghost.Eper, 'fro'), norm(ghost.Ebare, 'fro'), norm(ghost.Eimg, 'fro'), ...
                norm(ghost.parts.real, 'fro'), norm(ghost.parts.recip, 'fro'), norm(ghost.parts.surf, 'fro'), ...
                max(vecnorm(ghost.Eimg,2,2)), dE_eV, ...
                E_stationary_eV + dE_eV}; %#ok<SAGROW>

            fprintf(['    rcut %6.2f | kcut %5.2f | ||Eimg|| %.8e | ' ...
                     'max|Eimg| %.8e | dE %+0.8f eV | Ecorr %+0.8f eV\n'], ...
                rcutGhost, kcutGhost, norm(ghost.Eimg, 'fro'), ...
                max(vecnorm(ghost.Eimg,2,2)), dE_eV, E_stationary_eV + dE_eV);
        end
    end

    sweepRK = cell2table(rows, 'VariableNames', { ...
        'centerMode','alpha','rcut','kcut', ...
        'normEper','normEbare','normEimg', ...
        'normEreal','normErecip','normEsurf', ...
        'maxEimg','dEghostEV','EghostCorrectedEV'});

    fprintf('\n  Rcut/Kcut sweep table:\n');
    disp(sweepRK);

    %% --------------------------------------------------------------------
    % Alpha sweep at selected large ghost rcut/kcut

    fprintf('\n  Alpha sweep with ghost rcut %.2f and kcut %.2f:\n', ...
        cfg.ghostAlphaSweepRcut, cfg.ghostAlphaSweepKcut);

    nA = numel(cfg.ghostAlphaList);
    alphaRows = {};
    EimgAll = cell(nA, 1);

    for ia = 1:nA
        alphaGhost = cfg.ghostAlphaList(ia);

        ghost = local_point_dipole_image_field_real_geometry( ...
            targetPos, center, Mmu, Hrow, alphaGhost, ...
            cfg.ghostAlphaSweepRcut, cfg.ghostAlphaSweepKcut, cfg.scfBoundary);

        EimgAll{ia} = ghost.Eimg;

        dE_Ha = +0.5 * sum(sum(muPol .* ghost.Eimg));
        dE_eV = dE_Ha * 27.211386245988;

        alphaRows(end+1,:) = { ...
            centerMode, alphaGhost, cfg.ghostAlphaSweepRcut, cfg.ghostAlphaSweepKcut, ...
            norm(ghost.Eper, 'fro'), norm(ghost.Ebare, 'fro'), norm(ghost.Eimg, 'fro'), ...
            norm(ghost.parts.real, 'fro'), norm(ghost.parts.recip, 'fro'), norm(ghost.parts.surf, 'fro'), ...
            max(vecnorm(ghost.Eimg,2,2)), dE_eV, E_stationary_eV + dE_eV, NaN}; %#ok<SAGROW>

        fprintf(['    alpha %.3f | ||Ereal|| %.8e | ||Erecip|| %.8e | ' ...
                 '||Eimg|| %.8e | dE %+0.8f eV | Ecorr %+0.8f eV\n'], ...
            alphaGhost, norm(ghost.parts.real,'fro'), norm(ghost.parts.recip,'fro'), ...
            norm(ghost.Eimg,'fro'), dE_eV, E_stationary_eV + dE_eV);
    end

    EimgRef = EimgAll{end};
    for ia = 1:nA
        alphaRows{ia, end} = norm(EimgAll{ia} - EimgRef, 'fro') / ...
            max(1e-300, norm(EimgRef, 'fro'));
    end

    sweepAlpha = cell2table(alphaRows, 'VariableNames', { ...
        'centerMode','alpha','rcut','kcut', ...
        'normEper','normEbare','normEimg', ...
        'normEreal','normErecip','normEsurf', ...
        'maxEimg','dEghostEV','EghostCorrectedEV','relEimgToLast'});

    fprintf('\n  Alpha sweep table:\n');
    disp(sweepAlpha);

    %% --------------------------------------------------------------------
    % Print worst targets for final alpha setting

    ghostFinal = local_point_dipole_image_field_real_geometry( ...
        targetPos, center, Mmu, Hrow, cfg.ghostAlphaList(end), ...
        cfg.ghostAlphaSweepRcut, cfg.ghostAlphaSweepKcut, cfg.scfBoundary);

    EimgMag = vecnorm(ghostFinal.Eimg, 2, 2);
    EbareMag = vecnorm(ghostFinal.Ebare, 2, 2);
    EperMag = vecnorm(ghostFinal.Eper, 2, 2);

    [EimgSorted, EimgOrder] = sort(EimgMag, 'descend');

    fprintf('\n  Worst targets for final ghost settings:\n');
    fprintf('    alpha %.3f | rcut %.2f | kcut %.2f\n', ...
        cfg.ghostAlphaList(end), cfg.ghostAlphaSweepRcut, cfg.ghostAlphaSweepKcut);

    for kk = 1:min(cfg.nTopTargets, numel(EimgSorted))
        loc = EimgOrder(kk);
        site = targetSites(loc);

        fprintf(['    site %6d | dist_min %.8e | dist_raw %.8e | |Eimg| %.8e | |Eper| %.8e | ' ...
                 '|Ebare| %.8e | muDotEimg %.8e\n'], ...
            site, dMinMag(loc), dRawMag(loc), EimgMag(loc), EperMag(loc), EbareMag(loc), ...
            dot(muPol(loc,:), ghostFinal.Eimg(loc,:)));
    end
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

function center = local_choose_center(sys, sourceMask, Hrow, mode)
    switch lower(mode)
        case 'charged_mean'
            src = find(sourceMask);

            % Raw charged-site mean. If the charged molecules cross a
            % periodic boundary, this may still be a poor center convention.
            % The printed raw/min-image distances are meant to diagnose that.
            center = mean(sys.site_pos(src, :), 1);

        case 'cell_center'
            % Hrow stores direct lattice vectors as rows. The supercell
            % center from origin is half the sum of the row vectors.
            center = 0.5 * sum(Hrow, 1);

        otherwise
            error('Unknown center mode: %s', mode);
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
end

function ghost = local_point_dipole_image_field_real_geometry( ...
    targetPos, center, M, Hrow, alpha, rcut, kcut, boundary)

    % For a periodic point-dipole lattice, positions should be represented
    % modulo the supercell. The image-only subtraction
    %
    %   Eimg = Eper - Ebare
    %
    % must use the same minimum-image displacement for both Eper and Ebare.
    %
    % Hrow stores direct lattice vectors as rows:
    %   cart = frac * Hrow

    Hcol = Hrow.';

    dMin = local_minimum_image_displacements(targetPos, center, Hrow);

    targetLocal = dMin;
    centerLocal = [0.0, 0.0, 0.0];

    [Eper, parts] = local_periodic_point_dipole_field( ...
        targetLocal, centerLocal, M, Hcol, alpha, rcut, kcut, boundary);

    Ebare = local_bare_point_dipole_field(targetLocal, centerLocal, M);

    Eimg = Eper - Ebare;

    ghost = struct();
    ghost.center = center;
    ghost.M = M;
    ghost.alpha = alpha;
    ghost.rcut = rcut;
    ghost.kcut = kcut;
    ghost.boundary = boundary;
    ghost.dMin = dMin;
    ghost.targetLocal = targetLocal;
    ghost.Eper = Eper;
    ghost.Ebare = Ebare;
    ghost.Eimg = Eimg;
    ghost.parts = parts;
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
            error('Bad boundary.');
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