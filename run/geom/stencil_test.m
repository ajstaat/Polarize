%% benchmark_nonperiodic_directional_bin_size_sweep
% Sweep nonperiodic bin size on the active-space geometry used by the real
% nonperiodic workflow, and time the row-cache-like geometric workload.
%
% This script compares candidate bin sizes for the NONPERIODIC cell list.
% It keeps the geometry fixed and varies only opts.bin_size in
% geom.build_spatial_index(..., isPeriodic = false).
%
% Output focuses on:
%   - grid_shape
%   - search_reach
%   - n_offsets
%   - candidate pair count
%   - pair-set equality
%   - end-to-end query timing
%
% Notes:
%   - This benchmarks the nonperiodic spatial-index/query layer itself.
%   - It does not require the MEX row-cache builder to compare bin choices.
%   - If you want, a second version can benchmark build_active_row_cache
%     directly once we decide the best bin-size regime.

clear; clc; close all;

fprintf('=== nonperiodic directional bin-size sweep ===\n');

%% ------------------------------------------------------------------------
% Same real workflow as nonperiodic_real_crystal_example

cfg = struct();
cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename   = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

cfg.supercellSize = [3 7 2];
cfg.bondScale     = 1.20;
cfg.pairCharges   = [+1 -1];

cfg.use_thole  = true;
cfg.verbose    = true;
cfg.softening  = 0.0;
cfg.rcut_bohr  = 18.0;

cfg.nRepeat = 1;

% Sweep factors: bin_size = alpha * rcut
alphaList = [0.2, 0.225, 0.25, 0.275, 0.3, 0.325];

%% ------------------------------------------------------------------------
% Import crystal template

fprintf('\n[setup] importing crystal template...\n');
crystal = io.import_contcar_as_crystal(cfg.filename, ...
    'BondScale', cfg.bondScale, ...
    'WrapFractional', true, ...
    'SortMolecules', false);

fprintf('Imported crystal template:\n');
fprintf('  nSites    = %d\n', size(crystal.cart_coords, 1));
fprintf('  nBaseMols = %d\n', numel(unique(crystal.base_mol_id)));

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

fprintf('\nChosen pair (automatic):\n');
fprintf('  ref ID = %d\n', refID);
fprintf('  nbr ID = %d\n', nbrID);

%% ------------------------------------------------------------------------
% Apply charges and extract polsys

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
params_extract.ewald.mode = 'periodic';  % keep same extraction path as example

polsys = builder.extract_polarization_system(sys, params_extract);
io.assert_atomic_units(polsys);

fprintf('\nCanonical polarization system:\n');
fprintf('  n_sites             = %d\n', polsys.n_sites);
fprintf('  n_polarizable_sites = %d\n', nnz(polsys.site_is_polarizable));

%% ------------------------------------------------------------------------
% Build external field + active-space problem, same as real workflow

targetMask = logical(polsys.site_is_polarizable(:));
sourceMask = abs(polsys.site_charge(:)) > 0;

params_nonper = struct();
params_nonper.use_thole = cfg.use_thole;
params_nonper.field = struct();
params_nonper.field.mode = 'nonperiodic';
params_nonper.field.exclude_self = true;
params_nonper.field.use_thole_damping = cfg.use_thole;
params_nonper.field.target_mask = targetMask;
params_nonper.field.source_mask = sourceMask;

fprintf('\n[setup] building external field...\n');
Eext = calc.compute_external_field(polsys, params_nonper);

sorParams = struct();
sorParams.use_thole = cfg.use_thole;
sorParams.softening = cfg.softening;
sorParams.verbose = true;
sorParams.printEvery = 10;
sorParams.residualEvery = 10;
sorParams.tol = 1e-6;
sorParams.maxIter = 500;
sorParams.omega = 0.97;
sorParams.stopMetric = 'max_dmu';
sorParams.rcut = cfg.rcut_bohr;
sorParams.checkResidualAgainstLegacy = false;

fprintf('[setup] preparing SCF problem...\n');
problem = thole.prepare_scf_problem(polsys, Eext, sorParams);

posAct = polsys.site_pos(problem.activeSites, :);

fprintf('\nActive-space geometry:\n');
fprintf('  nActive = %d\n', size(posAct, 1));
fprintf('  rcut    = %.6f bohr\n', cfg.rcut_bohr);

%% ------------------------------------------------------------------------
% Reference case: current auto nonperiodic cell-list

optsRef = struct();
optsRef.isPeriodic = false;
optsRef.method     = 'cell_list';
optsRef.cutoff     = cfg.rcut_bohr;

tRefBase = tic;
spatial_ref = geom.build_spatial_index(posAct, optsRef);
timeRefBase = toc(tRefBase);

qopts = struct();
qopts.return_r        = true;
qopts.return_dr       = true;
qopts.return_full_idx = true;

pairs_ref = geom.query_pairs_within_cutoff(spatial_ref, cfg.rcut_bohr, qopts);
key_ref   = local_encode_pairs(pairs_ref.i, pairs_ref.j);

fprintf('\nReference (auto nonperiodic grid):\n');
fprintf('  grid_shape   = [%d %d %d]\n', spatial_ref.grid_shape);
fprintf('  search_reach = [%d %d %d]\n', spatial_ref.search_reach);
fprintf('  n_offsets    = %d\n', size(spatial_ref.neighbor_offsets, 1));
fprintf('  n_pairs      = %d\n', numel(pairs_ref.i));
fprintf('  build time   = %.6f s\n', timeRefBase);

%% ------------------------------------------------------------------------
% Sweep bin sizes

nCases = numel(alphaList);

gridShapeAll = zeros(nCases,3);
reachAll     = zeros(nCases,3);
nOffsetsAll  = zeros(nCases,1);
nPairsAll    = zeros(nCases,1);
buildAll     = zeros(nCases,1);
queryAll     = zeros(nCases,1);
matchPairs   = false(nCases,1);

for c = 1:nCases
    alpha = alphaList(c);
    binSize = alpha * cfg.rcut_bohr;

    opts = struct();
    opts.isPeriodic = false;
    opts.method     = 'cell_list';
    opts.cutoff     = cfg.rcut_bohr;
    opts.bin_size   = binSize;

    t0 = tic;
    spatial = geom.build_spatial_index(posAct, opts);
    buildAll(c) = toc(t0);

    pairs = geom.query_pairs_within_cutoff(spatial, cfg.rcut_bohr, qopts);
    key = local_encode_pairs(pairs.i, pairs.j);

    % Warmup
    geom.query_pairs_within_cutoff(spatial, cfg.rcut_bohr, qopts);

    tRep = zeros(cfg.nRepeat,1);
    for k = 1:cfg.nRepeat
        t1 = tic;
        geom.query_pairs_within_cutoff(spatial, cfg.rcut_bohr, qopts);
        tRep(k) = toc(t1);
    end

    gridShapeAll(c,:) = spatial.grid_shape;
    reachAll(c,:)     = spatial.search_reach;
    nOffsetsAll(c)    = size(spatial.neighbor_offsets, 1);
    nPairsAll(c)      = numel(pairs.i);
    queryAll(c)       = mean(tRep);
    matchPairs(c)     = isequal(key_ref, key);

    fprintf('\nalpha = %.2f  (bin_size = %.6f)\n', alpha, binSize);
    fprintf('  grid_shape   = [%d %d %d]\n', spatial.grid_shape);
    fprintf('  search_reach = [%d %d %d]\n', spatial.search_reach);
    fprintf('  n_offsets    = %d\n', nOffsetsAll(c));
    fprintf('  n_pairs      = %d\n', nPairsAll(c));
    fprintf('  pair-set match = %d\n', matchPairs(c));
    fprintf('  build        = %.6f s\n', buildAll(c));
    fprintf('  query        = %.6f s\n', queryAll(c));
end

%% ------------------------------------------------------------------------
% Summary

fprintf('\n============================================================\n');
fprintf('Summary\n');
fprintf('============================================================\n');
fprintf(' alpha  bin_size   grid_shape        reach      n_off  n_pairs  match  total(s)\n');

totalAll = buildAll + queryAll;
for c = 1:nCases
    fprintf(' %5.2f  %8.4f   [%3d %3d %3d]   [%d %d %d]   %5d  %7d    %d    %8.4f\n', ...
        alphaList(c), alphaList(c)*cfg.rcut_bohr, ...
        gridShapeAll(c,1), gridShapeAll(c,2), gridShapeAll(c,3), ...
        reachAll(c,1), reachAll(c,2), reachAll(c,3), ...
        nOffsetsAll(c), nPairsAll(c), matchPairs(c), totalAll(c));
end

validIdx = find(matchPairs);
if ~isempty(validIdx)
    [~, ii] = min(totalAll(validIdx));
    idxBest = validIdx(ii);

    fprintf('\nBest valid one-shot choice:\n');
    fprintf('  alpha        = %.2f\n', alphaList(idxBest));
    fprintf('  bin_size     = %.6f\n', alphaList(idxBest) * cfg.rcut_bohr);
    fprintf('  total time   = %.6f s\n', totalAll(idxBest));
    fprintf('  grid_shape   = [%d %d %d]\n', gridShapeAll(idxBest,:));
    fprintf('  search_reach = [%d %d %d]\n', reachAll(idxBest,:));
    fprintf('  n_offsets    = %d\n', nOffsetsAll(idxBest));
else
    fprintf('\nNo valid bin-size choices matched the reference result.\n');
end

fprintf('\nDone.\n');

%% ========================================================================
% Local helpers
%% ========================================================================

function key = local_encode_pairs(i, j)
key = uint64(i(:)) * uint64(2^32) + uint64(j(:));
key = sort(key);
end