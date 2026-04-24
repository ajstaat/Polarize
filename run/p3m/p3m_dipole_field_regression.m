%% run_p3m_dipole_field_regression
% Compare P3M reciprocal dipole field against direct reciprocal Ewald.
%
% First dipole milestone:
%
%   E_mu^recip,P3M ~= E_mu^recip,Ewald
%
% for a fixed dipole vector mu.
%
% This test does NOT compare real-space Thole damping, self terms, or total
% SCF energies. It isolates the long-range reciprocal dipole apply.

clear; clc; close all;

fprintf('=== P3M dipole reciprocal-field regression ===\n');

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

cfg.ewald = struct();
cfg.ewald.alpha = 0.35;
cfg.ewald.rcut = 11.5;
cfg.ewald.kcut = 3.5;
cfg.ewald.boundary = 'tinfoil';

% Mesh sweep.
meshList = { ...
    [24 40 24], ...
    [32 48 32], ...
    [40 64 40]};

orderList = [2 3 4];

% Fixed test dipole construction.
%
% 'linear_response_to_periodic_eext':
%   mu_i = alpha_i Eext_i using P3M periodic charge Eext.
%
% 'random':
%   deterministic random small dipoles on polarizable sites.
cfg.muMode = 'linear_response_to_periodic_eext';
cfg.randomMuScale = 0.03;  % au, used only for random mode

% P3M charge Eext settings used only to generate a realistic test mu.
cfg.eextP3M = struct();
cfg.eextP3M.mesh_size = [32 48 32];
cfg.eextP3M.assignment_order = 4;
cfg.eextP3M.realspace_backend = 'thole_periodic_real';
cfg.eextP3M.deconvolve_assignment = true;
cfg.eextP3M.deconvolution_floor = 1e-6;
cfg.eextP3M.derivative_mode = 'spectral';
cfg.eextP3M.influence_mode = 'ewald';
cfg.eextP3M.verbose = false;

%% ------------------------------------------------------------------------
% Import / build

fprintf('\n[setup] importing crystal template...\n');

crystal = io.import_contcar_as_crystal(cfg.filename, ...
    'BondScale', cfg.bondScale, ...
    'WrapFractional', true, ...
    'SortMolecules', false);

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

fprintf('\n[setup] building crystal system...\n');

buildOpts = struct();
buildOpts.supercell_size = cfg.supercellSize;
buildOpts.bondScale      = cfg.bondScale;
buildOpts.verbose        = cfg.verbose;
buildOpts.removeMolIDs   = [];
buildOpts.activeMolIDs   = [];

sys0 = builder.make_crystal_system(crystal, model, buildOpts);

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

fprintf('\n[setup] applying charges and extracting polsys...\n');

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
sourceChargeMask = abs(polsys.site_charge(:)) > 0;

fprintf('\nSystem summary:\n');
fprintf('  n_sites        = %d\n', polsys.n_sites);
fprintf('  n_targets      = %d\n', nnz(targetMask));
fprintf('  n_charge_sites = %d\n', nnz(sourceChargeMask));
fprintf('  total q source = %+0.16e\n', sum(polsys.site_charge(sourceChargeMask)));

%% ------------------------------------------------------------------------
% Build fixed test dipoles

fprintf('\n[setup] building fixed test dipole vector: %s\n', cfg.muMode);

mu = zeros(polsys.n_sites, 3);

switch lower(cfg.muMode)
    case 'linear_response_to_periodic_eext'
        p3mOpts = struct();
        p3mOpts.ewald = cfg.ewald;
        p3mOpts.mesh_size = cfg.eextP3M.mesh_size;
        p3mOpts.assignment_order = cfg.eextP3M.assignment_order;
        p3mOpts.target_mask = targetMask;
        p3mOpts.source_mask = sourceChargeMask;
        p3mOpts.exclude_self = true;
        p3mOpts.use_thole_damping = cfg.use_thole;
        p3mOpts.realspace_backend = cfg.eextP3M.realspace_backend;
        p3mOpts.deconvolve_assignment = cfg.eextP3M.deconvolve_assignment;
        p3mOpts.deconvolution_floor = cfg.eextP3M.deconvolution_floor;
        p3mOpts.derivative_mode = cfg.eextP3M.derivative_mode;
        p3mOpts.influence_mode = cfg.eextP3M.influence_mode;
        p3mOpts.verbose = cfg.eextP3M.verbose;

        [Eext, ~] = p3m.compute_external_field_charges(polsys, p3mOpts);

        alphaSite = polsys.site_alpha(:);
        mu(targetMask, :) = alphaSite(targetMask) .* Eext(targetMask, :);

    case 'random'
        rng(1);
        mu(targetMask, :) = cfg.randomMuScale .* randn(nnz(targetMask), 3);

    otherwise
        error('Unknown cfg.muMode: %s', cfg.muMode);
end

dipoleSourceMask = targetMask & (vecnorm(mu, 2, 2) > 0);

muPol = mu(targetMask, :);
muMag = vecnorm(muPol, 2, 2);

fprintf('  n_dipole_sources = %d\n', nnz(dipoleSourceMask));
fprintf('  ||mu||_F         = %.16e\n', norm(mu, 'fro'));
fprintf('  Mmu              = [%+.8e %+.8e %+.8e]\n', sum(mu(dipoleSourceMask,:), 1));
fprintf('  RMS |mu_i|       = %.16e au = %.8f D\n', ...
    sqrt(mean(muMag.^2)), sqrt(mean(muMag.^2))*2.541746473);
fprintf('  max |mu_i|       = %.16e au = %.8f D\n', ...
    max(muMag), max(muMag)*2.541746473);

%% ------------------------------------------------------------------------
% Direct reciprocal Ewald reference

fprintf('\n[reference] building direct reciprocal Ewald dipole field...\n');

tRef = tic;
[Eref, refParts] = local_direct_recip_dipole_field( ...
    polsys, mu, targetMask, dipoleSourceMask, cfg.ewald);
timeRef = toc(tRef);

ErefPol = Eref(targetMask, :);

Uref_Ha = -0.5 * sum(sum(muPol .* ErefPol));
Uref_eV = Uref_Ha * 27.211386245988;

fprintf('Reference done:\n');
fprintf('  time          = %.6f s\n', timeRef);
fprintf('  nK            = %d\n', refParts.nK);
fprintf('  ||Eref||_F    = %.16e\n', norm(Eref, 'fro'));
fprintf('  Urecip_ref    = %+0.16e Ha = %+0.8f eV\n', Uref_Ha, Uref_eV);

%% ------------------------------------------------------------------------
% P3M sweeps

rows = {};

fprintf('\n============================================================\n');
fprintf('P3M DIPOLE RECIPROCAL SWEEP\n');
fprintf('============================================================\n');

for im = 1:numel(meshList)
    meshSize = meshList{im};

    for io = 1:numel(orderList)
        order = orderList(io);

        fprintf('\n------------------------------------------------------------\n');
        fprintf('mesh = [%d %d %d], order = %d\n', meshSize, order);
        fprintf('------------------------------------------------------------\n');

        p3mOpts = struct();
        p3mOpts.ewald = cfg.ewald;
        p3mOpts.mesh_size = meshSize;
        p3mOpts.assignment_order = order;
        p3mOpts.target_mask = targetMask;
        p3mOpts.source_mask = dipoleSourceMask;
        p3mOpts.deconvolve_assignment = true;
        p3mOpts.deconvolution_floor = 1e-6;
        p3mOpts.verbose = true;

        tP3M = tic;
        [Ep3m, partsP3M] = p3m.compute_induced_field_dipoles( ...
            polsys, mu, p3mOpts);
        timeP3M = toc(tP3M);

        dEfield = Ep3m - Eref;

        targetSites = find(targetMask);
        errMag = vecnorm(dEfield(targetSites, :), 2, 2);
        refMag = vecnorm(Eref(targetSites, :), 2, 2);

        normRef = norm(Eref, 'fro');
        normP3M = norm(Ep3m, 'fro');
        normErr = norm(dEfield, 'fro');
        relErr = normErr / max(1e-300, normRef);

        rmsErr = sqrt(mean(errMag.^2));
        maxErr = max(errMag);
        rmsRef = sqrt(mean(refMag.^2));
        maxRef = max(refMag);

        Ep3mPol = Ep3m(targetMask, :);
        Up3m_Ha = -0.5 * sum(sum(muPol .* Ep3mPol));
        Up3m_eV = Up3m_Ha * 27.211386245988;

        dU_eV = Up3m_eV - Uref_eV;

        fprintf('\nP3M dipole reciprocal result:\n');
        fprintf('  time          = %.6f s\n', timeP3M);
        fprintf('  ||Ep3m||_F    = %.16e\n', normP3M);
        fprintf('  ||Eref||_F    = %.16e\n', normRef);
        fprintf('  ||diff||_F    = %.16e\n', normErr);
        fprintf('  rel diff      = %.16e\n', relErr);
        fprintf('  target RMS err= %.16e\n', rmsErr);
        fprintf('  target max err= %.16e\n', maxErr);
        fprintf('  Urecip_p3m    = %+0.16e Ha = %+0.8f eV\n', Up3m_Ha, Up3m_eV);
        fprintf('  dU_p3m-ref    = %+0.8e eV\n', dU_eV);

        rows(end+1,:) = { ...
            meshSize(1), meshSize(2), meshSize(3), order, ...
            normRef, normP3M, normErr, relErr, ...
            rmsRef, rmsErr, maxRef, maxErr, ...
            Uref_eV, Up3m_eV, dU_eV, ...
            partsP3M.nKmesh, ...
            timeRef, timeP3M, partsP3M.time_mesh}; %#ok<SAGROW>
    end
end

summary = cell2table(rows, 'VariableNames', { ...
    'M1','M2','M3','order', ...
    'normRef','normP3M','normErr','relErr', ...
    'rmsRef','rmsErr','maxRef','maxErr', ...
    'UrefEV','Up3mEV','dUeV', ...
    'nKmesh', ...
    'timeRef','timeP3M','timeMeshP3M'});

fprintf('\n============================================================\n');
fprintf('P3M DIPOLE REGRESSION SUMMARY\n');
fprintf('============================================================\n');
disp(summary);

fprintf('\nDone.\n');

%% ========================================================================
% Local helpers
%% ========================================================================

function [Erecip, parts] = local_direct_recip_dipole_field(sys, mu, targetMask, sourceMask, ew)
%LOCAL_DIRECT_RECIP_DIPOLE_FIELD Direct half-k reciprocal dipole Ewald field.
%
% This is a reference for the reciprocal part only:
%
%   E_rec(r_i) =
%      sum_{half k} two_pref(k)
%      [sum_j (k.mu_j) cos(k.(r_i-r_j))] k
%
% with
%
%   two_pref(k) = -8*pi/V exp(-k^2/(4 alpha^2)) / k^2
%
% This matches the reciprocal dipole tensor:
%
%   T_rec(k) = -4*pi/V exp(-k^2/(4 alpha^2)) k k^T / k^2
%
% using one representative from each +/-k pair.

    Hrow = local_get_direct_lattice(sys);
    Hcol = Hrow.';
    V = abs(det(Hcol));

    alpha = ew.alpha;
    kcut = ew.kcut;

    pos = sys.site_pos;

    targetSites = find(targetMask);
    sourceSites = find(sourceMask);

    targetPos = pos(targetSites, :);
    sourcePos = pos(sourceSites, :);
    muSrc = mu(sourceSites, :);

    Etarget = zeros(numel(targetSites), 3);

    [kvecs, meta] = ewald.enumerate_kvecs_triclinic(Hcol, kcut);
    nk = size(kvecs, 1);

    if nk > 0
        k2 = meta.k2(:);

        two_pref = -(8*pi/V) * exp(-k2 ./ (4 * alpha^2)) ./ k2;

        phase_src = sourcePos * kvecs.';
        cos_src = cos(phase_src);
        sin_src = sin(phase_src);

        mdotk = muSrc * kvecs.';  % nSources x nk

        rho_cos = sum(mdotk .* cos_src, 1);  % 1 x nk
        rho_sin = sum(mdotk .* sin_src, 1);  % 1 x nk

        phase_tgt = targetPos * kvecs.';
        cos_tgt = cos(phase_tgt);
        sin_tgt = sin(phase_tgt);

        amp = cos_tgt .* rho_cos + sin_tgt .* rho_sin;
        w = amp .* two_pref.';

        Etarget(:,1) = w * kvecs(:,1);
        Etarget(:,2) = w * kvecs(:,2);
        Etarget(:,3) = w * kvecs(:,3);
    end

    Erecip = zeros(size(sys.site_pos));
    Erecip(targetSites, :) = Etarget;

    parts = struct();
    parts.nK = nk;
    parts.alpha = alpha;
    parts.kcut = kcut;
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