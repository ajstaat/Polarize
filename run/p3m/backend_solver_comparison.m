%% backend_solver_comparison
% Compare:
%
%   A. Existing periodic SOR using:
%        - real-space row cache
%        - SOR/kCache reciprocal operator
%        - Ewald self/surface terms
%
%   B. Matrix-free GMRES/Krylov backend using:
%        - same real-space row cache
%        - cached P3M reciprocal dipole operator
%        - Ewald self/surface terms
%
% Requires:
%   p3m.build_dipole_cache
%   p3m.apply_dipole_cache
%
% Project lattice convention:
%
%   H rows are direct lattice vectors:
%
%       r_cart = f_frac * H
%
%   G columns are reciprocal lattice vectors:
%
%       H * G = 2*pi*I

clear; clc; close all;

fprintf('=== backend solver comparison: SOR/kCache vs cached P3M reciprocal GMRES ===\n');

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

% P3M fixed-charge external field.
cfg.p3mEext = struct();
cfg.p3mEext.mesh_size = [32 48 32];
cfg.p3mEext.assignment_order = 4;
cfg.p3mEext.realspace_backend = 'thole_periodic_real';
cfg.p3mEext.deconvolve_assignment = true;
cfg.p3mEext.deconvolution_floor = 1e-6;
cfg.p3mEext.derivative_mode = 'spectral';
cfg.p3mEext.influence_mode = 'ewald';
cfg.p3mEext.verbose = true;

% Cached P3M reciprocal dipole field used in GMRES matvec.
cfg.p3mDipole = struct();
cfg.p3mDipole.mesh_size = [40 64 40];
cfg.p3mDipole.assignment_order = 4;
cfg.p3mDipole.deconvolve_assignment = true;
cfg.p3mDipole.deconvolution_floor = 1e-6;
cfg.p3mDipole.use_kcut_mask = true;
cfg.p3mDipole.verbose = true;

% Existing periodic SOR reference.
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

% Matrix-free Krylov controls.
cfg.krylov = struct();
cfg.krylov.tol = 1e-6;
cfg.krylov.restart = 50;
cfg.krylov.maxOuter = 50;
cfg.krylov.verbose = true;

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

lat = geom.get_lattice(polsys);
lat2 = geom.get_lattice(lat);

fprintf('\nLattice convention:\n');
fprintf('  max |H*G - 2*pi*I|       = %.16e\n', ...
    max(abs(lat.H * lat.G - 2*pi*eye(3)), [], 'all'));
fprintf('  max |lat2.H*lat2.G - 2*pi*I| = %.16e\n', ...
    max(abs(lat2.H * lat2.G - 2*pi*eye(3)), [], 'all'));

targetMask = logical(polsys.site_is_polarizable(:));
sourceChargeMask = abs(polsys.site_charge(:)) > 0;

fprintf('\nSystem summary:\n');
fprintf('  n_sites        = %d\n', polsys.n_sites);
fprintf('  n_targets      = %d\n', nnz(targetMask));
fprintf('  n_charge_sites = %d\n', nnz(sourceChargeMask));
fprintf('  total q source = %+0.16e\n', sum(polsys.site_charge(sourceChargeMask)));

if nnz(targetMask & sourceChargeMask) ~= 0
    error('Charged source sites are still polarizable.');
end

%% ------------------------------------------------------------------------
% Build P3M periodic external field from fixed charges

fprintf('\n[Eext] building P3M periodic charge external field...\n');

p3mEextOpts = struct();
p3mEextOpts.ewald = cfg.ewald;
p3mEextOpts.mesh_size = cfg.p3mEext.mesh_size;
p3mEextOpts.assignment_order = cfg.p3mEext.assignment_order;
p3mEextOpts.target_mask = targetMask;
p3mEextOpts.source_mask = sourceChargeMask;
p3mEextOpts.exclude_self = true;
p3mEextOpts.use_thole_damping = cfg.use_thole;
p3mEextOpts.realspace_backend = cfg.p3mEext.realspace_backend;
p3mEextOpts.deconvolve_assignment = cfg.p3mEext.deconvolve_assignment;
p3mEextOpts.deconvolution_floor = cfg.p3mEext.deconvolution_floor;
p3mEextOpts.derivative_mode = cfg.p3mEext.derivative_mode;
p3mEextOpts.influence_mode = cfg.p3mEext.influence_mode;
p3mEextOpts.verbose = cfg.p3mEext.verbose;

tEext = tic;
[Eext, eextParts] = p3m.compute_external_field_charges(polsys, p3mEextOpts);
timeEext = toc(tEext);

fprintf('\nEext summary:\n');
fprintf('  time        = %.6f s\n', timeEext);
fprintf('  ||Eext||_F  = %.16e\n', norm(Eext, 'fro'));
fprintf('  ||Ereal||_F = %.16e\n', norm(eextParts.real, 'fro'));
fprintf('  ||Erecip||_F= %.16e\n', norm(eextParts.recip, 'fro'));
fprintf('  ||Esurf||_F = %.16e\n', norm(eextParts.surf, 'fro'));
if isfield(eextParts, 'rho_total')
    fprintf('  sum rho mesh= %+0.16e\n', eextParts.rho_total);
end

%% ------------------------------------------------------------------------
% Prepare SCF problem, real-space row cache, and kCache

fprintf('\n[setup] preparing SCF problem...\n');

sorParams = cfg.sor;
sorParams.use_thole = cfg.use_thole;
sorParams.kspace_mode = cfg.sor.kspace_mode;
sorParams.k_block_size = cfg.sor.k_block_size;
sorParams.kspace_memory_limit_gb = cfg.sor.kspace_memory_limit_gb;
sorParams.use_mex_kspace = cfg.sor.use_mex_kspace;

problem = thole.prepare_scf_problem(polsys, Eext, sorParams);

fprintf('  nPolSites = %d\n', problem.nPolSites);

fprintf('\n[setup] building periodic real-space dipole row cache...\n');

rowOpts = struct();
rowOpts.profile = true;
rowOpts.use_mex = cfg.sor.use_mex_realspace;
rowOpts.use_thole = cfg.use_thole;

tRow = tic;
rowCache = geom.build_active_row_cache_periodic( ...
    polsys, problem, cfg.ewald, rowOpts);
timeRow = toc(tRow);

fprintf('  row cache time  = %.6f s\n', timeRow);
fprintf('  row cache count = %d\n', local_row_cache_count(rowCache));

fprintf('\n[setup] building SOR k-space cache...\n');

kOpts = struct();
kOpts.kspace_mode = cfg.sor.kspace_mode;
kOpts.kspace_memory_limit_gb = cfg.sor.kspace_memory_limit_gb;
kOpts.k_block_size = cfg.sor.k_block_size;
kOpts.verbose = true;

tK = tic;
kCache = geom.build_periodic_kspace_cache(polsys, problem, cfg.ewald, kOpts);
timeK = toc(tK);

fprintf('  k cache time = %.6f s\n', timeK);
fprintf('  nK           = %d\n', kCache.num_kvec);
fprintf('  storage mode = %s\n', kCache.storage_mode);
if isfield(kCache, 'lattice_convention')
    fprintf('  convention   = %s\n', kCache.lattice_convention);
end
if isfield(kCache, 'G')
    fprintf('  max |H*G - 2*pi*I| = %.16e\n', ...
        max(abs(kCache.H * kCache.G - 2*pi*eye(3)), [], 'all'));
end

sorParams.realspace_dipole_row_cache = rowCache;
sorParams.kspace_cache = kCache;

%% ------------------------------------------------------------------------
% Build P3M dipole cache

fprintf('\n[setup] building P3M dipole cache...\n');

p3mDipoleCacheOpts = struct();
p3mDipoleCacheOpts.ewald = cfg.ewald;
p3mDipoleCacheOpts.mesh_size = cfg.p3mDipole.mesh_size;
p3mDipoleCacheOpts.assignment_order = cfg.p3mDipole.assignment_order;
p3mDipoleCacheOpts.target_mask = targetMask;
p3mDipoleCacheOpts.source_mask = targetMask;
p3mDipoleCacheOpts.deconvolve_assignment = cfg.p3mDipole.deconvolve_assignment;
p3mDipoleCacheOpts.deconvolution_floor = cfg.p3mDipole.deconvolution_floor;
p3mDipoleCacheOpts.use_kcut_mask = cfg.p3mDipole.use_kcut_mask;
p3mDipoleCacheOpts.verbose = cfg.p3mDipole.verbose;

tP3MCache = tic;
p3mDipoleCache = p3m.build_dipole_cache(polsys, p3mDipoleCacheOpts);
timeP3MCache = toc(tP3MCache);

fprintf('  P3M dipole cache time = %.6f s\n', timeP3MCache);
fprintf('  P3M dipole cache nK   = %d\n', p3mDipoleCache.nK);

%% ------------------------------------------------------------------------
% Reference existing periodic SOR

fprintf('\n============================================================\n');
fprintf('REFERENCE: existing periodic SOR\n');
fprintf('============================================================\n');

tRef = tic;
[muRef, infoRef] = thole.solve_scf_iterative_periodic_sor( ...
    polsys, Eext, cfg.ewald, sorParams);
timeRef = toc(tRef);

[ErefHa, ErefEV] = local_stationary_energy(targetMask, muRef, Eext);

fprintf('\nReference SOR result:\n');
fprintf('  time      = %.6f s\n', timeRef);
fprintf('  converged = %d\n', infoRef.converged);
fprintf('  nIter     = %d\n', infoRef.nIter);
fprintf('  relres    = %.16e\n', infoRef.relres);
fprintf('  ||mu||_F  = %.16e\n', norm(muRef, 'fro'));
fprintf('  energy    = %+0.16e Ha = %+0.8f eV\n', ErefHa, ErefEV);

%% ------------------------------------------------------------------------
% Matvec diagnostic at reference SOR mu

fprintf('\n============================================================\n');
fprintf('MATVEC DIAGNOSTIC AT REFERENCE SOR MU\n');
fprintf('============================================================\n');

diagRef = local_operator_diagnostic_at_mu( ...
    polsys, lat, Eext, targetMask, rowCache, p3mDipoleCache, muRef, cfg.ewald);

fprintf('Operator field norms at muRef:\n');
fprintf('  ||EopNeeded||_F       = %.16e\n', diagRef.normNeeded);
fprintf('  ||ErealLocal||_F      = %.16e\n', diagRef.normReal);
fprintf('  ||ErecipP3M||_F       = %.16e\n', diagRef.normRecip);
fprintf('  ||Eself||_F           = %.16e\n', diagRef.normSelf);
fprintf('  ||Esurf||_F           = %.16e\n', diagRef.normSurf);
fprintf('  ||EopLocal||_F        = %.16e\n', diagRef.normLocal);
fprintf('  ||EopLocal-needed||_F = %.16e\n', diagRef.normDiff);
fprintf('  rel operator error    = %.16e\n', diagRef.relErr);

fprintf('\nSign checks:\n');
fprintf('  dot(needed, local) / norms = %.16e\n', diagRef.cosLocal);
fprintf('  dot(needed, real)  / norms = %.16e\n', diagRef.cosReal);
fprintf('  dot(needed, recip) / norms = %.16e\n', diagRef.cosRecip);
fprintf('  dot(needed, self)  / norms = %.16e\n', diagRef.cosSelf);
fprintf('  dot(needed, surf)  / norms = %.16e\n', diagRef.cosSurf);

fprintf('\nEnergy contractions at muRef:\n');
fprintf('  -1/2 mu.EopNeeded = %+0.8f eV\n', diagRef.UneededEV);
fprintf('  -1/2 mu.Ereal     = %+0.8f eV\n', diagRef.UrealEV);
fprintf('  -1/2 mu.Erecip    = %+0.8f eV\n', diagRef.UrecipEV);
fprintf('  -1/2 mu.Eself     = %+0.8f eV\n', diagRef.UselfEV);
fprintf('  -1/2 mu.Esurf     = %+0.8f eV\n', diagRef.UsurfEV);
fprintf('  -1/2 mu.EopLocal  = %+0.8f eV\n', diagRef.UlocalEV);

%% ------------------------------------------------------------------------
% P3M reciprocal backend Krylov solve

fprintf('\n============================================================\n');
fprintf('Cached P3M reciprocal backend Krylov solve\n');
fprintf('============================================================\n');

tP3M = tic;
[muP3M, infoP3M] = local_solve_scf_periodic_p3m_gmres( ...
    polsys, lat, Eext, problem, rowCache, p3mDipoleCache, cfg.ewald, cfg.krylov);
timeP3M = toc(tP3M);

[Ep3mHa, Ep3mEV] = local_stationary_energy(targetMask, muP3M, Eext);

fprintf('\nCached P3M Krylov result:\n');
fprintf('  time      = %.6f s\n', timeP3M);
fprintf('  converged = %d\n', infoP3M.converged);
fprintf('  flag      = %d\n', infoP3M.flag);
fprintf('  iter      = %s\n', mat2str(infoP3M.iter));
fprintf('  relres    = %.16e\n', infoP3M.relres);
fprintf('  fpRelres  = %.16e\n', infoP3M.fixedPointRelres);
fprintf('  matvecs   = %d\n', infoP3M.nMatvec);
fprintf('  ||mu||_F  = %.16e\n', norm(muP3M, 'fro'));
fprintf('  energy    = %+0.16e Ha = %+0.8f eV\n', Ep3mHa, Ep3mEV);

%% ------------------------------------------------------------------------
% Comparison

muDiff = muP3M - muRef;

polSites = find(targetMask);
muRefPol = muRef(polSites, :);
muDiffPol = muDiff(polSites, :);

normMuRef = norm(muRef, 'fro');
normMuDiff = norm(muDiff, 'fro');
relMuDiff = normMuDiff / max(1e-300, normMuRef);

siteDiffMag = vecnorm(muDiffPol, 2, 2);
siteRefMag = vecnorm(muRefPol, 2, 2);

fprintf('\n============================================================\n');
fprintf('BACKEND COMPARISON SUMMARY\n');
fprintf('============================================================\n');

fprintf('External field / setup:\n');
fprintf('  P3M Eext time                 = %.6f s\n', timeEext);
fprintf('  row cache build time          = %.6f s\n', timeRow);
fprintf('  kCache build time             = %.6f s\n', timeK);
fprintf('  P3M dipole cache build time   = %.6f s\n', timeP3MCache);

fprintf('\nReference SOR:\n');
fprintf('  time                          = %.6f s\n', timeRef);
fprintf('  converged / nIter / relres    = %d / %d / %.3e\n', ...
    infoRef.converged, infoRef.nIter, infoRef.relres);
fprintf('  energy                        = %+0.8f eV\n', ErefEV);

fprintf('\nCached P3M Krylov:\n');
fprintf('  time                          = %.6f s\n', timeP3M);
fprintf('  converged / iter / relres     = %d / %s / %.3e\n', ...
    infoP3M.converged, mat2str(infoP3M.iter), infoP3M.relres);
fprintf('  fixed-point relres            = %.3e\n', infoP3M.fixedPointRelres);
fprintf('  matvecs                       = %d\n', infoP3M.nMatvec);
fprintf('  energy                        = %+0.8f eV\n', Ep3mEV);

fprintf('\nDifferences:\n');
fprintf('  dE(P3M - ref)                 = %+0.8e eV\n', Ep3mEV - ErefEV);
fprintf('  ||muP3M - muRef||_F           = %.16e\n', normMuDiff);
fprintf('  rel ||dmu||                   = %.16e\n', relMuDiff);
fprintf('  RMS site |dmu|                = %.16e au = %.8e D\n', ...
    sqrt(mean(siteDiffMag.^2)), sqrt(mean(siteDiffMag.^2))*2.541746473);
fprintf('  max site |dmu|                = %.16e au = %.8e D\n', ...
    max(siteDiffMag), max(siteDiffMag)*2.541746473);
fprintf('  RMS site |muRef|              = %.16e au = %.8e D\n', ...
    sqrt(mean(siteRefMag.^2)), sqrt(mean(siteRefMag.^2))*2.541746473);

fprintf('\nCached P3M timing breakdown:\n');
fprintf('  total real apply time         = %.6f s\n', infoP3M.timeRealApply);
fprintf('  total P3M reciprocal time     = %.6f s\n', infoP3M.timeP3MRecip);
fprintf('  total self/surface time       = %.6f s\n', infoP3M.timeSelfSurf);
fprintf('  avg real apply / matvec       = %.6f s\n', ...
    infoP3M.timeRealApply / max(1, infoP3M.nMatvec));
fprintf('  avg P3M reciprocal / matvec   = %.6f s\n', ...
    infoP3M.timeP3MRecip / max(1, infoP3M.nMatvec));

fprintf('\nDone.\n');

%% ========================================================================
% Local solver/helper functions
%% ========================================================================

function [mu, info] = local_solve_scf_periodic_p3m_gmres( ...
    sys, lat, Eext, problem, rowCache, p3mDipoleCache, ewaldParams, krylovOpts)

    polSites = local_get_problem_pol_sites(problem, sys);
    nPol = numel(polSites);

    alphaPol = sys.site_alpha(polSites);
    EextPol = Eext(polSites, :);

    bPol = alphaPol .* EextPol;
    b = bPol(:);

    alphaEwald = ewaldParams.alpha;
    selfCoeff = +(4 * alphaEwald^3 / (3 * sqrt(pi)));

    surfCoeff = local_dipole_surface_coeff(lat, ewaldParams);

    tol = local_get_opt(krylovOpts, 'tol', 1e-6);
    restart = local_get_opt(krylovOpts, 'restart', 50);
    maxOuter = local_get_opt(krylovOpts, 'maxOuter', 50);
    verbose = local_get_opt(krylovOpts, 'verbose', true);

    nMatvec = 0;
    timeRealApply = 0.0;
    timeP3MRecip = 0.0;
    timeSelfSurf = 0.0;

    if verbose
        fprintf('GMRES(cached P3M periodic backend): nPol=%d, nDOF=%d, restart=%d, maxOuter=%d, tol=%.3e\n', ...
            nPol, 3*nPol, restart, maxOuter, tol);
        fprintf('  Ewald dipole selfCoeff = %.16e\n', selfCoeff);
        fprintf('  Ewald dipole surfCoeff = %.16e\n', surfCoeff);
        fprintf('  P3M cache mesh         = [%d %d %d]\n', p3mDipoleCache.mesh_size);
        fprintf('  P3M cache order        = %d\n', p3mDipoleCache.assignment_order);
        fprintf('  P3M cache nK           = %d\n', p3mDipoleCache.nK);
    end

    Afun = @(x) local_matvec_A(x);

    tSolve = tic;
    [x, flag, relres, iter, resvec] = gmres(Afun, b, restart, tol, maxOuter);
    timeSolve = toc(tSolve);

    muPol = reshape(x, [nPol, 3]);

    mu = zeros(sys.n_sites, 3);
    mu(polSites, :) = muPol;

    [~, EopPolFinal] = local_matvec_A(x);
    EtotalPol = EextPol + EopPolFinal;
    fpResid = muPol - alphaPol .* EtotalPol;
    fixedPointRelres = norm(fpResid, 'fro') / max(1e-300, norm(muPol, 'fro'));

    info = struct();
    info.converged = (flag == 0);
    info.flag = flag;
    info.relres = relres;
    info.iter = iter;
    info.resvec = resvec;
    info.fixedPointRelres = fixedPointRelres;
    info.nMatvec = nMatvec;
    info.timeSolve = timeSolve;
    info.timeRealApply = timeRealApply;
    info.timeP3MRecip = timeP3MRecip;
    info.timeSelfSurf = timeSelfSurf;
    info.avgRealApply = timeRealApply / max(1, nMatvec);
    info.avgP3MRecip = timeP3MRecip / max(1, nMatvec);
    info.selfCoeff = selfCoeff;
    info.surfCoeff = surfCoeff;

    if verbose
        fprintf('GMRES done:\n');
        fprintf('  flag       = %d\n', flag);
        fprintf('  relres     = %.16e\n', relres);
        fprintf('  iter       = %s\n', mat2str(iter));
        fprintf('  matvecs    = %d\n', nMatvec);
        fprintf('  time solve = %.6f s\n', timeSolve);
        fprintf('  time real  = %.6f s total | %.6f s/matvec\n', ...
            timeRealApply, info.avgRealApply);
        fprintf('  time p3m   = %.6f s total | %.6f s/matvec\n', ...
            timeP3MRecip, info.avgP3MRecip);
        fprintf('  fpRelres   = %.16e\n', fixedPointRelres);
    end

    function [y, EopPol] = local_matvec_A(xvec)
        nMatvec = nMatvec + 1;

        muTrialPol = reshape(xvec, [nPol, 3]);

        tReal = tic;
        ErealPol = local_apply_dipole_row_cache(rowCache, muTrialPol);
        timeRealApply = timeRealApply + toc(tReal);

        muTrialFull = zeros(sys.n_sites, 3);
        muTrialFull(polSites, :) = muTrialPol;

        tP3M = tic;
        [ErecipFull, ~] = p3m.apply_dipole_cache(p3mDipoleCache, muTrialFull);
        timeP3MRecip = timeP3MRecip + toc(tP3M);

        ErecipPol = ErecipFull(polSites, :);

        tSS = tic;
        EselfPol = selfCoeff .* muTrialPol;

        if surfCoeff ~= 0
            Mmu = sum(muTrialPol, 1);
            EsurfPol = repmat(surfCoeff .* Mmu, nPol, 1);
        else
            EsurfPol = zeros(nPol, 3);
        end
        timeSelfSurf = timeSelfSurf + toc(tSS);

        EopPol = ErealPol + ErecipPol + EselfPol + EsurfPol;

        yPol = muTrialPol - alphaPol .* EopPol;
        y = yPol(:);
    end
end

function diag = local_operator_diagnostic_at_mu( ...
    sys, lat, Eext, targetMask, rowCache, p3mDipoleCache, mu, ewaldParams)

    polSites = find(targetMask);
    alphaPol = sys.site_alpha(polSites);

    muPol = mu(polSites, :);
    EextPol = Eext(polSites, :);

    alphaEwald = ewaldParams.alpha;
    selfCoeff = +(4 * alphaEwald^3 / (3 * sqrt(pi)));
    surfCoeff = local_dipole_surface_coeff(lat, ewaldParams);

    EopNeededPol = muPol ./ alphaPol - EextPol;

    ErealPol = local_apply_dipole_row_cache(rowCache, muPol);

    [ErecipFull, ~] = p3m.apply_dipole_cache(p3mDipoleCache, mu);
    ErecipPol = ErecipFull(polSites, :);

    EselfPol = selfCoeff .* muPol;

    if surfCoeff ~= 0
        Mmu = sum(muPol, 1);
        EsurfPol = repmat(surfCoeff .* Mmu, numel(polSites), 1);
    else
        EsurfPol = zeros(numel(polSites), 3);
    end

    EopLocalPol = ErealPol + ErecipPol + EselfPol + EsurfPol;

    dOp = EopLocalPol - EopNeededPol;

    diag = struct();
    diag.EopNeededPol = EopNeededPol;
    diag.ErealPol = ErealPol;
    diag.ErecipPol = ErecipPol;
    diag.EselfPol = EselfPol;
    diag.EsurfPol = EsurfPol;
    diag.EopLocalPol = EopLocalPol;
    diag.dOp = dOp;

    diag.normNeeded = norm(EopNeededPol, 'fro');
    diag.normReal = norm(ErealPol, 'fro');
    diag.normRecip = norm(ErecipPol, 'fro');
    diag.normSelf = norm(EselfPol, 'fro');
    diag.normSurf = norm(EsurfPol, 'fro');
    diag.normLocal = norm(EopLocalPol, 'fro');
    diag.normDiff = norm(dOp, 'fro');
    diag.relErr = diag.normDiff / max(1e-300, diag.normNeeded);

    diag.cosLocal = local_cosine(EopNeededPol, EopLocalPol);
    diag.cosReal = local_cosine(EopNeededPol, ErealPol);
    diag.cosRecip = local_cosine(EopNeededPol, ErecipPol);
    diag.cosSelf = local_cosine(EopNeededPol, EselfPol);
    diag.cosSurf = local_cosine(EopNeededPol, EsurfPol);

    Ha_to_eV = 27.211386245988;
    diag.UneededEV = -0.5 * sum(sum(muPol .* EopNeededPol)) * Ha_to_eV;
    diag.UrealEV = -0.5 * sum(sum(muPol .* ErealPol)) * Ha_to_eV;
    diag.UrecipEV = -0.5 * sum(sum(muPol .* ErecipPol)) * Ha_to_eV;
    diag.UselfEV = -0.5 * sum(sum(muPol .* EselfPol)) * Ha_to_eV;
    diag.UsurfEV = -0.5 * sum(sum(muPol .* EsurfPol)) * Ha_to_eV;
    diag.UlocalEV = -0.5 * sum(sum(muPol .* EopLocalPol)) * Ha_to_eV;

    diag.selfCoeff = selfCoeff;
    diag.surfCoeff = surfCoeff;
end

function Epol = local_apply_dipole_row_cache(rowCache, muPol)

    row_ptr = rowCache.row_ptr(:);
    col_idx = rowCache.col_idx(:);

    nRows = numel(row_ptr) - 1;
    Epol = zeros(nRows, 3);

    if ~(isfield(rowCache, 'coeff_iso') && ...
         isfield(rowCache, 'coeff_dyad') && ...
         isfield(rowCache, 'dr'))
        fprintf('\nrowCache fields are:\n');
        disp(fieldnames(rowCache));
        error('Expected rowCache.coeff_iso, coeff_dyad, and dr.');
    end

    coeff_iso_all = rowCache.coeff_iso(:);
    coeff_dyad_all = rowCache.coeff_dyad(:);
    dr_all = rowCache.dr;

    for i = 1:nRows
        idx0 = row_ptr(i);
        idx1 = row_ptr(i + 1) - 1;

        if idx1 < idx0
            continue;
        end

        idx = idx0:idx1;
        cols = col_idx(idx);

        muNbr = muPol(cols, :);
        dr = dr_all(idx, :);

        muDotR = sum(muNbr .* dr, 2);
        contrib = coeff_iso_all(idx) .* muNbr + ...
                  coeff_dyad_all(idx) .* (muDotR .* dr);

        Epol(i, :) = sum(contrib, 1);
    end
end

function surfCoeff = local_dipole_surface_coeff(lat, ewaldParams)

    boundary = 'tinfoil';
    if isfield(ewaldParams, 'boundary') && ~isempty(ewaldParams.boundary)
        boundary = lower(char(string(ewaldParams.boundary)));
    end

    switch boundary
        case 'tinfoil'
            surfCoeff = 0.0;

        case 'vacuum'
            surfCoeff = -4*pi/(3*lat.volume);

        otherwise
            error('Unknown boundary: %s', boundary);
    end
end

function polSites = local_get_problem_pol_sites(problem, sys)
    if isfield(problem, 'activeSites') && ~isempty(problem.activeSites)
        polSites = problem.activeSites(:);
        return;
    end

    candidates = {'polSites','pol_sites','active_sites','site_indices','pol_idx','active_idx'};

    for k = 1:numel(candidates)
        name = candidates{k};
        if isfield(problem, name) && ~isempty(problem.(name))
            polSites = problem.(name)(:);
            return;
        end
    end

    polSites = find(logical(sys.site_is_polarizable(:)));
end

function [EHa, EeV] = local_stationary_energy(targetMask, mu, Eext)
    polSites = find(targetMask);
    muPol = mu(polSites, :);
    EextPol = Eext(polSites, :);

    EHa = -0.5 * sum(sum(muPol .* EextPol));
    EeV = EHa * 27.211386245988;
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

function c = local_cosine(A, B)
    denom = norm(A, 'fro') * norm(B, 'fro');
    if denom < 1e-300
        c = NaN;
    else
        c = sum(A(:) .* B(:)) / denom;
    end
end

function val = local_get_opt(s, name, defaultVal)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = defaultVal;
    end
end