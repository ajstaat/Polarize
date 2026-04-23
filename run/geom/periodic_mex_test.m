%% real_periodic_rowcache_validation_medium
% Validate periodic real-space row-cache on a medium real crystal by comparing:
%   1) brute-force MIC reference
%   2) MATLAB periodic row-cache path
%   3) MEX periodic row-cache path
%
% This test activates all polarizable sites and compares the real-space row cache.
% It also reconstructs UNIQUE undirected pair sets from each row cache so we can
% diagnose whether mismatches are due to:
%   - duplicate insertion, or
%   - genuinely wrong neighbor ownership / missing pairs.

clear; clc; close all;

fprintf('=== real periodic row-cache validation (medium system) ===\n');

%% ------------------------------------------------------------------------
% User controls

cfg = struct();

cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename   = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

cfg.supercellSize = [3 6 2];
cfg.bondScale     = 1.20;
cfg.verbose       = true;

cfg.periodic = struct();
cfg.periodic.alpha    = 0.30;
cfg.periodic.rcut     = 11.5;
cfg.periodic.kcut     = 3.5;
cfg.periodic.boundary = 'tinfoil';

cfg.use_thole = true;

cfg.compareTolAbs = 1e-11;
cfg.compareTolRel = 1e-9;

cfg.maxExamplePairs = 20;

ewaldParams = struct();
ewaldParams.alpha    = cfg.periodic.alpha;
ewaldParams.rcut     = cfg.periodic.rcut;
ewaldParams.kcut     = cfg.periodic.kcut;
ewaldParams.boundary = cfg.periodic.boundary;

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

sys = builder.make_crystal_system(crystal, model, buildOpts);

fprintf('\nBuilt working system:\n');
fprintf('  n_unit_sites = %d\n', sys.n_unit_sites);
fprintf('  n_cells      = %d\n', sys.n_cells);
fprintf('  n_sites      = %d\n', sys.n_sites);

%% ------------------------------------------------------------------------
% Extract canonical polarization system

fprintf('\n[setup] extracting canonical polarization system...\n');
params_extract = struct();
params_extract.ewald = struct();
params_extract.ewald.mode = 'periodic';

polsys = builder.extract_polarization_system(sys, params_extract);
io.assert_atomic_units(polsys);

fprintf('Canonical polarization system:\n');
fprintf('  n_sites             = %d\n', polsys.n_sites);
fprintf('  n_polarizable_sites = %d\n', nnz(polsys.site_is_polarizable));

%% ------------------------------------------------------------------------
% Problem = all polarizable sites active

problem = struct();
problem.activeSites = find(polsys.site_is_polarizable(:));

fprintf('\nProblem summary:\n');
fprintf('  nActive = %d\n', numel(problem.activeSites));

%% ------------------------------------------------------------------------
% Lattice diagnostics

H = local_get_direct_lattice(polsys);
Lmin = geom.shortest_lattice_translation(H);
rcutMaxSafe = 0.5 * Lmin;

fprintf('\nPeriodic real-space geometry:\n');
fprintf('  lattice lengths = [%.6f %.6f %.6f] bohr\n', ...
    norm(H(:,1)), norm(H(:,2)), norm(H(:,3)));
fprintf('  Lmin           = %.6f bohr\n', Lmin);
fprintf('  Lmin/2         = %.6f bohr\n', rcutMaxSafe);
fprintf('  requested rcut = %.6f bohr\n', cfg.periodic.rcut);

if cfg.periodic.rcut >= rcutMaxSafe
    error(['Chosen rcut violates the single-image condition for this cell.\n' ...
           'Require rcut < Lmin/2 = %.6f bohr, but got rcut = %.6f bohr.'], ...
           rcutMaxSafe, cfg.periodic.rcut);
end

%% ------------------------------------------------------------------------
% Build periodic cell-list spatial index explicitly

fprintf('\n[setup] building periodic cell-list spatial index explicitly...\n');
posAct = polsys.site_pos(problem.activeSites, :);

spatialOpts = struct();
spatialOpts.isPeriodic = true;
spatialOpts.cell = H;
spatialOpts.method = 'cell_list';
spatialOpts.cutoff = cfg.periodic.rcut;

tSpatial = tic;
spatial = geom.build_spatial_index(posAct, spatialOpts);
timeSpatial = toc(tSpatial);

fprintf('[setup] spatial backend = %s\n', spatial.backend);
fprintf('        build time = %.6f s\n', timeSpatial);
fprintf('        grid_shape = [%d %d %d]\n', spatial.grid_shape);
fprintf('        search_reach = [%d %d %d]\n', spatial.search_reach);
fprintf('        n_offsets = %d\n', size(spatial.neighbor_offsets,1));

% %% ------------------------------------------------------------------------
% % Brute-force reference
% 
% fprintf('\n[ref] building brute-force periodic MIC reference...\n');
% tRef = tic;
% ref = local_build_reference_rowcache(polsys, problem, ewaldParams, cfg.use_thole);
% timeRef = toc(tRef);
% 
% fprintf('[ref] done in %.6f s\n', timeRef);
% fprintf('  nPairsUndirected = %d\n', ref.nPairsUndirected);
% fprintf('  nEntriesDirected = %d\n', ref.nEntriesDirected);

% %% ------------------------------------------------------------------------
% % MATLAB path
% 
% fprintf('\n[matlab] building row cache...\n');
% optsMat = struct();
% optsMat.profile = true;
% optsMat.use_mex = false;
% optsMat.use_thole = cfg.use_thole;
% 
% tMat = tic;
% rowMat = geom.build_active_row_cache_periodic(polsys, problem, ewaldParams, optsMat, spatial);
% timeMat = toc(tMat);
% 
% fprintf('[matlab] done in %.6f s\n', timeMat);
% fprintf('  nPairsUndirected = %d\n', rowMat.nPairsUndirected);
% fprintf('  nEntriesDirected = %d\n', rowMat.nEntriesDirected);

%% ------------------------------------------------------------------------
% MEX path

% mexExists = (exist(['mex_build_active_row_cache_periodic.' mexext], 'file') == 3) || ...
%             (exist('mex_build_active_row_cache_periodic', 'file') == 3);
% 
% if ~mexExists
%     error(['Periodic MEX backend was not found on the MATLAB path.\n' ...
%            'Compile mex_build_active_row_cache_periodic first, then rerun this test.']);
% end

fprintf('\n[mex] building row cache...\n');
optsMex = struct();
optsMex.profile = true;
optsMex.use_mex = true;
optsMex.use_thole = cfg.use_thole;

tMex = tic;
rowMex = geom.build_active_row_cache_periodic(polsys, problem, ewaldParams, optsMex, spatial);
timeMex = toc(tMex);

fprintf('[mex] done in %.6f s\n', timeMex);
fprintf('  nPairsUndirected = %d\n', rowMex.nPairsUndirected);
fprintf('  nEntriesDirected = %d\n', rowMex.nEntriesDirected);

%% ------------------------------------------------------------------------
% Timing summary

fprintf('\n============================================================\n');
fprintf('Timing summary\n');
fprintf('============================================================\n');
fprintf('  spatial index build : %.6f s\n', timeSpatial);
%fprintf('  brute-force ref     : %.6f s\n', timeRef);
%fprintf('  MATLAB row cache    : %.6f s\n', timeMat);
fprintf('  MEX row cache       : %.6f s\n', timeMex);
%fprintf('  speedup (MAT/MEX)   : %.3f x\n', timeMat / timeMex);

% %% ------------------------------------------------------------------------
% % Raw cache comparisons
% 
% fprintf('\n============================================================\n');
% fprintf('Reference vs MATLAB\n');
% fprintf('============================================================\n');
% local_compare_rowcaches(ref, rowMat, cfg.compareTolAbs, cfg.compareTolRel);
% 
% fprintf('\n============================================================\n');
% fprintf('Reference vs MEX\n');
% fprintf('============================================================\n');
% local_compare_rowcaches(ref, rowMex, cfg.compareTolAbs, cfg.compareTolRel);
% 
% %% ------------------------------------------------------------------------
% % Unique undirected pair-set diagnostics
% 
% fprintf('\n============================================================\n');
% fprintf('Unique undirected pair-set diagnostics\n');
% fprintf('============================================================\n');
% 
% Pref = geom.rowcache_undirected_pair_set(ref);
% Pmat = geom.rowcache_undirected_pair_set(rowMat);
% Pmex = geom.rowcache_undirected_pair_set(rowMex);
% 
% fprintf('  unique undirected pairs (REF) = %d\n', size(Pref,1));
% fprintf('  unique undirected pairs (MAT) = %d\n', size(Pmat,1));
% fprintf('  unique undirected pairs (MEX) = %d\n', size(Pmex,1));
% 
% fprintf('  duplicate directed entries implied:\n');
% fprintf('    REF : %d\n', ref.nPairsUndirected - size(Pref,1));
% fprintf('    MAT : %d\n', rowMat.nPairsUndirected - size(Pmat,1));
% fprintf('    MEX : %d\n', rowMex.nPairsUndirected - size(Pmex,1));
% 
% extraMat = setdiff(Pmat, Pref, 'rows');
% missMat  = setdiff(Pref, Pmat, 'rows');
% 
% extraMex = setdiff(Pmex, Pref, 'rows');
% missMex  = setdiff(Pref, Pmex, 'rows');
% 
% fprintf('\nMATLAB pair-set differences vs REF:\n');
% fprintf('  extra pairs = %d\n', size(extraMat,1));
% fprintf('  missing pairs = %d\n', size(missMat,1));
% 
% fprintf('\nMEX pair-set differences vs REF:\n');
% fprintf('  extra pairs = %d\n', size(extraMex,1));
% fprintf('  missing pairs = %d\n', size(missMex,1));
% 
% if ~isempty(extraMat)
%     fprintf('\nExample extra MATLAB pairs (first %d):\n', min(cfg.maxExamplePairs, size(extraMat,1)));
%     disp(extraMat(1:min(cfg.maxExamplePairs, size(extraMat,1)), :));
% end
% if ~isempty(missMat)
%     fprintf('\nExample missing MATLAB pairs (first %d):\n', min(cfg.maxExamplePairs, size(missMat,1)));
%     disp(missMat(1:min(cfg.maxExamplePairs, size(missMat,1)), :));
% end
% if ~isempty(extraMex)
%     fprintf('\nExample extra MEX pairs (first %d):\n', min(cfg.maxExamplePairs, size(extraMex,1)));
%     disp(extraMex(1:min(cfg.maxExamplePairs, size(extraMex,1)), :));
% end
% if ~isempty(missMex)
%     fprintf('\nExample missing MEX pairs (first %d):\n', min(cfg.maxExamplePairs, size(missMex,1)));
%     disp(missMex(1:min(cfg.maxExamplePairs, size(missMex,1)), :));
% end

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

function rowCache = local_build_reference_rowcache(sys, problem, ewaldParams, use_thole)
    pos = sys.site_pos;
    activeSites = problem.activeSites(:);
    posAct = pos(activeSites, :);
    nActive = numel(activeSites);

    H = local_get_direct_lattice(sys);
    Hinv = inv(H);

    alphaEwald = ewaldParams.alpha;
    rcut = ewaldParams.rcut;
    cutoff2 = rcut^2;

    pairI = [];
    pairJ = [];
    drUndir = zeros(0,3);
    r2Undir = zeros(0,1);
    rUndir  = zeros(0,1);
    coeffIsoUndir = zeros(0,1);
    coeffDyadUndir = zeros(0,1);

    for i = 1:nActive
        for j = i+1:nActive
            drRaw = posAct(j,:) - posAct(i,:);
            drMic = local_mic_displacement(drRaw, H, Hinv);
            r2 = sum(drMic.^2);

            if r2 <= cutoff2
                r = sqrt(r2);

                invR2 = 1 / r2;
                invR  = 1 / r;
                invR3 = invR * invR2;
                invR5 = invR3 * invR2;
                invR4 = invR2^2;

                alpha2 = alphaEwald^2;
                twoAlphaOverSqrtPi = 2 * alphaEwald / sqrt(pi);

                erfcar = erfc(alphaEwald * r);
                expar2 = exp(-alpha2 * r2);

                B = erfcar * invR3 + twoAlphaOverSqrtPi * expar2 * invR2;
                C = 3 * erfcar * invR5 + ...
                    twoAlphaOverSqrtPi * (2 * alpha2 * invR2 + 3 * invR4) * expar2;

                coeffIso = -B;
                coeffDyad = +C;

                if use_thole
                    ai = sys.site_alpha(activeSites(i));
                    aj = sys.site_alpha(activeSites(j));
                    [f3, f5] = local_thole_f3f5_scalar(r, ai, aj, sys.thole_a);
                    l3 = f3 - 1;
                    l5 = f5 - 1;
                    coeffIso = coeffIso - l3 * invR3;
                    coeffDyad = coeffDyad + 3 * l5 * invR5;
                end

                pairI(end+1,1) = i; %#ok<AGROW>
                pairJ(end+1,1) = j; %#ok<AGROW>
                drUndir(end+1,:) = drMic; %#ok<AGROW>
                r2Undir(end+1,1) = r2; %#ok<AGROW>
                rUndir(end+1,1)  = r; %#ok<AGROW>
                coeffIsoUndir(end+1,1) = coeffIso; %#ok<AGROW>
                coeffDyadUndir(end+1,1) = coeffDyad; %#ok<AGROW>
            end
        end
    end

    nPairs = numel(pairI);
    rowCounts = accumarray([pairI; pairJ], 1, [nActive,1], @sum, 0);

    row_ptr = zeros(nActive+1,1);
    row_ptr(1) = 1;
    row_ptr(2:end) = 1 + cumsum(rowCounts);

    nDir = 2*nPairs;
    col_idx = zeros(nDir,1);
    drDir = zeros(nDir,3);
    r2Dir = zeros(nDir,1);
    rDir  = zeros(nDir,1);
    coeffIsoDir = zeros(nDir,1);
    coeffDyadDir = zeros(nDir,1);

    nextPtr = row_ptr(1:end-1);

    for p = 1:nPairs
        i = pairI(p);
        j = pairJ(p);

        k = nextPtr(i);
        col_idx(k) = j;
        drDir(k,:) = drUndir(p,:);
        r2Dir(k) = r2Undir(p);
        rDir(k) = rUndir(p);
        coeffIsoDir(k) = coeffIsoUndir(p);
        coeffDyadDir(k) = coeffDyadUndir(p);
        nextPtr(i) = k + 1;

        k = nextPtr(j);
        col_idx(k) = i;
        drDir(k,:) = -drUndir(p,:);
        r2Dir(k) = r2Undir(p);
        rDir(k) = rUndir(p);
        coeffIsoDir(k) = coeffIsoUndir(p);
        coeffDyadDir(k) = coeffDyadUndir(p);
        nextPtr(j) = k + 1;
    end

    rowCache = struct();
    rowCache.activeSites = activeSites;
    rowCache.nActive = nActive;
    rowCache.row_ptr = row_ptr;
    rowCache.col_idx = col_idx;
    rowCache.dr = drDir;
    rowCache.r2_bare = r2Dir;
    rowCache.r_bare = rDir;
    rowCache.coeff_iso = coeffIsoDir;
    rowCache.coeff_dyad = coeffDyadDir;
    rowCache.nPairsUndirected = nPairs;
    rowCache.nEntriesDirected = nDir;
end

function drMic = local_mic_displacement(drRaw, H, Hinv)
    sf = (Hinv * drRaw(:)).';
    sf = sf - round(sf);
    drMic = (H * sf(:)).';
end

function [f3, f5] = local_thole_f3f5_scalar(r, alpha_i, alpha_j, a)
    if r <= 0 || alpha_i <= 0 || alpha_j <= 0 || a <= 0
        f3 = 1;
        f5 = 1;
        return;
    end
    aij = (alpha_i * alpha_j)^(1/6);
    s = a * (r / aij)^3;
    e = exp(-s);
    f3 = 1 - e;
    f5 = 1 - (1 + s) * e;
end

function local_compare_rowcaches(A, B, absTol, relTol)
    fprintf('  nPairsUndirected match = %d\n', A.nPairsUndirected == B.nPairsUndirected);
    fprintf('  nEntriesDirected match = %d\n', A.nEntriesDirected == B.nEntriesDirected);
    fprintf('  row_ptr exact match    = %d\n', isequal(A.row_ptr, B.row_ptr));
    fprintf('  col_idx exact match    = %d\n', isequal(A.col_idx, B.col_idx));

    local_report_array_compare('dr', A.dr, B.dr, absTol, relTol);
    local_report_array_compare('r2_bare', A.r2_bare, B.r2_bare, absTol, relTol);
    local_report_array_compare('r_bare', A.r_bare, B.r_bare, absTol, relTol);
    local_report_array_compare('coeff_iso', A.coeff_iso, B.coeff_iso, absTol, relTol);
    local_report_array_compare('coeff_dyad', A.coeff_dyad, B.coeff_dyad, absTol, relTol);
end

function local_report_array_compare(name, A, B, absTol, relTol)
    if ~isequal(size(A), size(B))
        fprintf('  %s size match      = 0\n', name);
        fprintf('    size(A) = [%s]\n', num2str(size(A)));
        fprintf('    size(B) = [%s]\n', num2str(size(B)));
        return;
    end

    d = abs(A(:) - B(:));
    maxAbs = max(d);

    denom = max(max(abs(A(:))), max(abs(B(:))));
    if denom == 0
        maxRel = 0;
    else
        maxRel = maxAbs / denom;
    end

    ok = (maxAbs <= absTol) || (maxRel <= relTol);

    fprintf('  %s within tol      = %d\n', name, ok);
    fprintf('    max abs diff     = %.6e\n', maxAbs);
    fprintf('    max rel diff     = %.6e\n', maxRel);
end