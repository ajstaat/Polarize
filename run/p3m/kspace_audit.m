%% run_lattice_reciprocal_convention_audit
% Audit lattice / reciprocal-lattice conventions after introducing:
%
%   lat = geom.get_lattice(sys)
%
% Project convention:
%
%   H has direct lattice vectors as rows:
%
%       r_cart = f_frac * H
%
%   G has reciprocal lattice vectors as columns:
%
%       H * G = 2*pi*I
%
%   For integer reciprocal index m_col:
%
%       k_col = G * m_col
%       k_row = k_col.'
%
% This script checks:
%
%   1. geom.get_lattice identity H*G = 2*pi*I
%   2. ewald.enumerate_kvecs_from_lattice matches kCache
%   3. p3m.make_kgrid uses the same project convention
%   4. legacy enumerate_kvecs_triclinic(Hrow) / (Hrow.') behavior is shown
%      only as a comparison, not as the source of truth.
%
% This script does NOT run SCF.

clear; clc; close all;

fprintf('=== lattice / reciprocal convention audit ===\n');

%% ------------------------------------------------------------------------
% User controls

cfg = struct();

cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename   = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

cfg.supercellSize = [2 5 1];
cfg.bondScale     = 1.20;
cfg.verbose       = true;

cfg.ewald = struct();
cfg.ewald.alpha = 0.35;
cfg.ewald.rcut  = 11.5;
cfg.ewald.kcut  = 3.5;
cfg.ewald.boundary = 'tinfoil';

cfg.meshSize = [32 48 32];

setTol = 1e-10;

%% ------------------------------------------------------------------------
% Import / build system

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

sys = builder.make_crystal_system(crystal, model, buildOpts);

params_extract = struct();
params_extract.ewald = struct();
params_extract.ewald.mode = 'periodic';

polsys = builder.extract_polarization_system(sys, params_extract);
io.assert_atomic_units(polsys);

lat = geom.get_lattice(polsys);

H = lat.H;
G = lat.G;
V = lat.volume;

fprintf('\n============================================================\n');
fprintf('PROJECT LATTICE CONVENTION\n');
fprintf('============================================================\n');

fprintf('\nProject convention:\n');
fprintf('  H rows are direct lattice vectors.\n');
fprintf('  G columns are reciprocal lattice vectors.\n');
fprintf('  r_cart = f_frac * H\n');
fprintf('  k_col  = G * m_col\n');
fprintf('  H * G  = 2*pi*I\n');

fprintf('\nH from geom.get_lattice(polsys):\n');
disp(H);

fprintf('G from geom.get_lattice(polsys):\n');
disp(G);

fprintf('det(H) / volume = %.16e / %.16e\n', det(H), V);

fprintf('\nDirect lattice vector lengths, norm(H(i,:)):\n');
disp(vecnorm(H, 2, 2).');

fprintf('Reciprocal lattice vector lengths, norm(G(:,i)):\n');
disp(vecnorm(G, 2, 1));

Id = H * G;

fprintf('\nH * G:\n');
disp(Id);

fprintf('max abs error from 2*pi*I = %.16e\n', ...
    max(abs(Id - 2*pi*eye(3)), [], 'all'));

assert(max(abs(Id - 2*pi*eye(3)), [], 'all') < 1e-10, ...
    'Project reciprocal identity H*G = 2*pi*I failed.');

%% ------------------------------------------------------------------------
% Legacy comparison only

fprintf('\n============================================================\n');
fprintf('LEGACY ENUMERATOR COMPARISON ONLY\n');
fprintf('============================================================\n');

Hrow = H;
Hcol = H.';

fprintf('\nLegacy column-style enumerator expects direct lattice vectors as columns.\n');
fprintf('These two calls are shown only to expose old/new convention differences.\n');

[k_legacy_Hrow, meta_legacy_Hrow] = ewald.enumerate_kvecs_triclinic(Hrow, cfg.ewald.kcut);
[k_legacy_Hcol, meta_legacy_Hcol] = ewald.enumerate_kvecs_triclinic(Hcol, cfg.ewald.kcut);

fprintf('\nlegacy enumerate_kvecs_triclinic(Hrow, kcut):\n');
fprintf('  nK       = %d\n', size(k_legacy_Hrow, 1));
fprintf('  hkmax    = [%d %d %d]\n', meta_legacy_Hrow.hkmax);
fprintf('  min |k|  = %.16e\n', min(meta_legacy_Hrow.knorm));
fprintf('  max |k|  = %.16e\n', max(meta_legacy_Hrow.knorm));

fprintf('\nlegacy enumerate_kvecs_triclinic(Hrow.'', kcut):\n');
fprintf('  nK       = %d\n', size(k_legacy_Hcol, 1));
fprintf('  hkmax    = [%d %d %d]\n', meta_legacy_Hcol.hkmax);
fprintf('  min |k|  = %.16e\n', min(meta_legacy_Hcol.knorm));
fprintf('  max |k|  = %.16e\n', max(meta_legacy_Hcol.knorm));

fprintf('\nCompare legacy Hrow vs legacy Hrow.'':\n');
cmp_legacy = local_compare_ksets(k_legacy_Hrow, k_legacy_Hcol, setTol);
local_print_kset_comparison(cmp_legacy);

%% ------------------------------------------------------------------------
% Project-convention enumeration

fprintf('\n============================================================\n');
fprintf('PROJECT ENUMERATOR CHECK\n');
fprintf('============================================================\n');

[k_project, meta_project] = ewald.enumerate_kvecs_from_lattice(lat, cfg.ewald.kcut);

fprintf('\newald.enumerate_kvecs_from_lattice(lat, kcut):\n');
fprintf('  nK       = %d\n', size(k_project, 1));
fprintf('  hkmax    = [%d %d %d]\n', meta_project.hkmax);
fprintf('  min |k|  = %.16e\n', min(meta_project.knorm));
fprintf('  max |k|  = %.16e\n', max(meta_project.knorm));

fprintf('\nFirst few project k-vectors:\n');
disp(k_project(1:min(10,size(k_project,1)), :));

fprintf('\nCompare project enumerator to legacy enumerate(Hrow):\n');
cmp_project_vs_legacy_Hrow = local_compare_ksets(k_project, k_legacy_Hrow, setTol);
local_print_kset_comparison(cmp_project_vs_legacy_Hrow);

fprintf('\nCompare project enumerator to legacy enumerate(Hrow.''):\n');
cmp_project_vs_legacy_Hcol = local_compare_ksets(k_project, k_legacy_Hcol, setTol);
local_print_kset_comparison(cmp_project_vs_legacy_Hcol);

%% ------------------------------------------------------------------------
% Production kCache check

fprintf('\n============================================================\n');
fprintf('PRODUCTION kCache CHECK\n');
fprintf('============================================================\n');

sorParams = struct();
sorParams.use_thole = true;
sorParams.verbose = false;
sorParams.tol = 1e-6;
sorParams.maxIter = 1;
sorParams.omega = 0.55;
sorParams.stopMetric = 'max_dmu';

Ezero = zeros(polsys.n_sites, 3);
problem = thole.prepare_scf_problem(polsys, Ezero, sorParams);

kOpts = struct();
kOpts.kspace_mode = 'full';
kOpts.kspace_memory_limit_gb = 8;
kOpts.k_block_size = 2048;
kOpts.verbose = true;

kCache = geom.build_periodic_kspace_cache(polsys, problem, cfg.ewald, kOpts);

fprintf('\nkCache summary:\n');
fprintf('  nK            = %d\n', kCache.num_kvec);
fprintf('  storage_mode  = %s\n', kCache.storage_mode);
fprintf('  kcut          = %.16e\n', kCache.kcut);
fprintf('  min |k|       = %.16e\n', min(kCache.knorm));
fprintf('  max |k|       = %.16e\n', max(kCache.knorm));

if isfield(kCache, 'lattice_convention')
    fprintf('  convention    = %s\n', kCache.lattice_convention);
else
    fprintf('  convention    = <missing>\n');
end

if isfield(kCache, 'G')
    fprintf('\nkCache.H * kCache.G max error from 2*pi*I = %.16e\n', ...
        max(abs(kCache.H * kCache.G - 2*pi*eye(3)), [], 'all'));
end

fprintf('\nCompare kCache.kvecs to project enumerator:\n');
cmp_cache_project = local_compare_ksets(kCache.kvecs, k_project, setTol);
local_print_kset_comparison(cmp_cache_project);

fprintf('\nCompare kCache.kvecs to legacy enumerate(Hrow):\n');
cmp_cache_legacy_Hrow = local_compare_ksets(kCache.kvecs, k_legacy_Hrow, setTol);
local_print_kset_comparison(cmp_cache_legacy_Hrow);

fprintf('\nCompare kCache.kvecs to legacy enumerate(Hrow.''):\n');
cmp_cache_legacy_Hcol = local_compare_ksets(kCache.kvecs, k_legacy_Hcol, setTol);
local_print_kset_comparison(cmp_cache_legacy_Hcol);

assert(cmp_cache_project.nCommon == cmp_cache_project.nA && ...
       cmp_cache_project.nCommon == cmp_cache_project.nB, ...
       'kCache.kvecs does not exactly match project enumerator.');

%% ------------------------------------------------------------------------
% P3M kgrid check

fprintf('\n============================================================\n');
fprintf('P3M KGRID CHECK\n');
fprintf('============================================================\n');

kg = p3m.make_kgrid(lat, cfg.meshSize);

mask = (kg.k2 > 0) & (kg.k2 <= cfg.ewald.kcut^2);
kgrid_full = [kg.kx(mask), kg.ky(mask), kg.kz(mask)];

% Compare the P3M full +/- set to the full +/- version of the project
% enumerator. This avoids arbitrary half-space selection in Cartesian
% components.
kproject_full = [k_project; -k_project];

fprintf('\nP3M kgrid using p3m.make_kgrid(lat, meshSize):\n');
fprintf('  meshSize           = [%d %d %d]\n', cfg.meshSize);
fprintf('  nK full +/-        = %d\n', size(kgrid_full, 1));
fprintf('  min |k| full       = %.16e\n', min(vecnorm(kgrid_full,2,2)));
fprintf('  max |k| full       = %.16e\n', max(vecnorm(kgrid_full,2,2)));

fprintf('\nCompare P3M full kgrid to full +/- project enumerator:\n');
cmp_p3m_project_full = local_compare_ksets(kgrid_full, kproject_full, setTol);
local_print_kset_comparison(cmp_p3m_project_full);

fprintf('\nInterpretation:\n');
fprintf('  nCommon < nProjectFull is expected for finite FFT mesh support.\n');
fprintf('  Zero overlap is a convention failure. Substantial overlap is good.\n');

assert(cmp_p3m_project_full.nCommon > 0, ...
    'P3M kgrid has zero overlap with project reciprocal convention.');

%% ------------------------------------------------------------------------
% Phase sanity check

fprintf('\n============================================================\n');
fprintf('PHASE SANITY CHECK\n');
fprintf('============================================================\n');

pos = polsys.site_pos;
activeSites = problem.activeSites(:);
posPol = pos(activeSites, :);

nTestSites = min(5, size(posPol, 1));
nTestK = min(5, size(k_project, 1));

fprintf('\nFor Cartesian pos and Cartesian k, phase = pos * k.''.\n');
fprintf('Sample phases using project enumerator first few k:\n');

phaseSample = posPol(1:nTestSites, :) * k_project(1:nTestK, :).';
disp(phaseSample);

fprintf('\nFractional phase sanity:\n');

fracPol = posPol(1:nTestSites, :) / H;
hkl = meta_project.hkl(1:nTestK, :);
phaseFrac = 2*pi * fracPol * hkl.';
phaseCart = phaseSample;

fprintf('max |phaseCart - phaseFrac| for sample = %.16e\n', ...
    max(abs(phaseCart - phaseFrac), [], 'all'));

disp('phaseCart:');
disp(phaseCart);

disp('phaseFrac:');
disp(phaseFrac);

assert(max(abs(phaseCart - phaseFrac), [], 'all') < 1e-8, ...
    'Cartesian phase and fractional phase disagree.');

fprintf('\n============================================================\n');
fprintf('AUDIT SUMMARY\n');
fprintf('============================================================\n');

fprintf('Project identity H*G passed.\n');
fprintf('kCache matches ewald.enumerate_kvecs_from_lattice(lat,kcut).\n');
fprintf('P3M kgrid has nonzero overlap with project reciprocal convention.\n');
fprintf('Cartesian and fractional phase agree for sampled modes.\n');

fprintf('\nDone.\n');

%% ========================================================================
% Local helpers
%% ========================================================================

function cmp = local_compare_ksets(A, B, tol)
% Compare two k-vector sets up to rounding tolerance.
%
% A and B are both Nx3 arrays of Cartesian row k-vectors.

    A = local_unique_kset(local_round_kset(A, tol), tol);
    B = local_unique_kset(local_round_kset(B, tol), tol);

    keyA = local_kset_keys(A);
    keyB = local_kset_keys(B);

    [commonKeys, ia, ib] = intersect(keyA, keyB, 'stable'); %#ok<ASGLU>
    [missingKeys, imiss] = setdiff(keyB, keyA, 'stable'); %#ok<ASGLU>
    [extraKeys, iextra] = setdiff(keyA, keyB, 'stable'); %#ok<ASGLU>

    cmp = struct();
    cmp.nA = size(A, 1);
    cmp.nB = size(B, 1);
    cmp.nCommon = numel(commonKeys);
    cmp.nMissingFromA = numel(missingKeys); % in B not A
    cmp.nExtraInA = numel(extraKeys);       % in A not B
    cmp.fracCommonOfA = cmp.nCommon / max(1, cmp.nA);
    cmp.fracCommonOfB = cmp.nCommon / max(1, cmp.nB);

    if ~isempty(commonKeys)
        Acommon = A(ia, :);
        Bcommon = B(ib, :);
        d = Acommon - Bcommon;
        cmp.maxAbsCommonDiff = max(abs(d), [], 'all');
        cmp.rmsCommonDiff = sqrt(mean(d(:).^2));
    else
        cmp.maxAbsCommonDiff = NaN;
        cmp.rmsCommonDiff = NaN;
    end

    cmp.sampleExtraInA = A(iextra(1:min(10,numel(iextra))), :);
    cmp.sampleMissingFromA = B(imiss(1:min(10,numel(imiss))), :);
end

function local_print_kset_comparison(cmp)
    fprintf('  nA                   = %d\n', cmp.nA);
    fprintf('  nB                   = %d\n', cmp.nB);
    fprintf('  nCommon              = %d\n', cmp.nCommon);
    fprintf('  nExtraInA            = %d\n', cmp.nExtraInA);
    fprintf('  nMissingFromA        = %d\n', cmp.nMissingFromA);
    fprintf('  fracCommonOfA        = %.8f\n', cmp.fracCommonOfA);
    fprintf('  fracCommonOfB        = %.8f\n', cmp.fracCommonOfB);
    fprintf('  maxAbsCommonDiff     = %.16e\n', cmp.maxAbsCommonDiff);
    fprintf('  rmsCommonDiff        = %.16e\n', cmp.rmsCommonDiff);

    if ~isempty(cmp.sampleExtraInA)
        fprintf('  sample extra in A:\n');
        disp(cmp.sampleExtraInA);
    end
    if ~isempty(cmp.sampleMissingFromA)
        fprintf('  sample missing from A:\n');
        disp(cmp.sampleMissingFromA);
    end
end

function K = local_round_kset(K, tol)
    K = round(K ./ tol) .* tol;
    K(abs(K) < tol/2) = 0;
end

function K = local_unique_kset(K, tol)
    K = local_round_kset(K, tol);
    keys = local_kset_keys(K);
    [~, ia] = unique(keys, 'stable');
    K = K(ia, :);
end

function keys = local_kset_keys(K)
    keys = strings(size(K,1), 1);
    for i = 1:size(K,1)
        keys(i) = sprintf('%+.12e,%+.12e,%+.12e', K(i,1), K(i,2), K(i,3));
    end
end