%% run_vasp_periodic_dielectric_tensor_workflow
%
% Compute induced-dipole / electronic dielectric tensor for a neutral
% periodic molecular crystal.
%
% Workflow:
%
%   VASP/CONTCAR
%   -> neutral periodic crystal system
%   -> periodic polarization system
%   -> build periodic P3M polarization operator once
%   -> apply uniform Cartesian fields Ex, Ey, Ez
%   -> solve induced dipoles
%   -> compute chi and epsilon:
%
%        P = M / V
%        chi(:,j) = P / E0_j
%        epsilon = I + 4*pi*chi
%
% Notes:
%
%   - This computes the clamped-geometry induced-dipole dielectric response
%     of the polarizable-site model, not an ionic/relaxed-ion dielectric.
%   - The applied field is a uniform Cartesian external field.
%   - The tensor is reported in Cartesian coordinates.
%   - For this linear model, the response should be independent of field
%     amplitude as long as the solver is converged.

clear; clc; close all;

fprintf('\n============================================================\n');
fprintf('Periodic induced-dipole dielectric tensor workflow\n');
fprintf('============================================================\n');

%% ------------------------------------------------------------------------
% User controls
% -------------------------------------------------------------------------

cfg = struct();

cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

% For dielectric tensor, use a neutral periodic cell. This does not need to
% be the same huge cell used for pair finite-size convergence, but it should
% be large enough that the molecular/template construction is robust.
cfg.supercellSize = [3 11 3];
cfg.bondScale = 1.20;

cfg.verbose = true;

% Uniform field amplitude in atomic units.
%
% 1 a.u. field = 5.14220674763e11 V/m.
% 1e-4 a.u. is already enormous experimentally, but this is a linear model,
% so amplitude only affects numerical scaling.
cfg.field_amplitude_au = 1e-5;

% Optional: also run -E for a finite-difference/symmetry check.
cfg.run_plus_minus_check = false;

% Periodic P3M settings.
cfg.ewald = struct();
cfg.ewald.alpha = 0.30;
cfg.ewald.rcut = 18.0;
cfg.ewald.kcut = 2.50;

% For dielectric response, tinfoil is the usual clean periodic bulk setting.
% You can test 'vacuum' separately, but it changes the macroscopic boundary
% condition rather than the local polarizable model.
cfg.ewald.boundary = 'tinfoil';

cfg.p3m = struct();

% Mesh suggestion for [3 11 3]. Adjust if you change supercell size.
cfg.p3m.mesh_size = [64 216 64];
cfg.p3m.assignment_order = 4;

cfg.p3m.derivative_mode = 'spectral';
cfg.p3m.influence_mode = 'ewald';
cfg.p3m.fd_stencil = 'central2';
cfg.p3m.deconvolve_assignment = true;
cfg.p3m.deconvolution_floor = 1e-8;
cfg.p3m.alias_range = 2;

% Solver settings.
cfg.scf = struct();
cfg.scf.tol = 1e-8;
cfg.scf.maxIter = 500;
cfg.scf.verbose = false;

cfg.solver = struct();
cfg.solver.method = 'gmres';
cfg.solver.gmres_restart = [];
cfg.solver.compute_residual = true;

% Operator settings.
cfg.operator = struct();
cfg.operator.backend = 'auto';
cfg.operator.use_thole = true;
cfg.operator.softening = 0.0;
cfg.operator.use_mex = true;
cfg.operator.profile = false;
cfg.operator.verbose = false;

%% ------------------------------------------------------------------------
% Check input
% -------------------------------------------------------------------------

if ~isfile(cfg.filename)
    error('VASP file not found:\n  %s\nEdit cfg.filename at the top of this script.', cfg.filename);
end

fprintf('\nInput structure:\n  %s\n', cfg.filename);

fprintf('\nRun controls:\n');
fprintf('  supercell              = [%d %d %d]\n', cfg.supercellSize);
fprintf('  field amplitude        = %.6e a.u.\n', cfg.field_amplitude_au);
fprintf('  plus/minus check       = %d\n', cfg.run_plus_minus_check);
fprintf('  alpha                  = %.6g\n', cfg.ewald.alpha);
fprintf('  rcut                   = %.6g bohr\n', cfg.ewald.rcut);
fprintf('  kcut                   = %.6g bohr^-1\n', cfg.ewald.kcut);
fprintf('  boundary               = %s\n', cfg.ewald.boundary);
fprintf('  P3M mesh               = [%d %d %d]\n', cfg.p3m.mesh_size);
fprintf('  P3M assignment order   = %d\n', cfg.p3m.assignment_order);
fprintf('  solver                 = %s\n', cfg.solver.method);
fprintf('  tolerance              = %.3e\n', cfg.scf.tol);

%% ------------------------------------------------------------------------
% 1. Import crystal template
% -------------------------------------------------------------------------

fprintf('\n[1] Importing crystal template...\n');

tImport = tic;

crystal = io.import_contcar_as_crystal(cfg.filename, ...
    'BondScale', cfg.bondScale, ...
    'SortMolecules', false);

importTime = toc(tImport);

fprintf('  import time = %.6f s\n', importTime);
fprintf('  nSites      = %d\n', crystal.nSites);
fprintf('  nBaseMols   = %d\n', crystal.nBaseMols);
fprintf('  units       = %s\n', crystal.units.length);

%% ------------------------------------------------------------------------
% 2. Define polarizability model
% -------------------------------------------------------------------------

fprintf('\n[2] Defining polarizability model...\n');

model = struct();
model.thole_a = 0.39;
model.alpha_units = 'angstrom^3';

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
    'C_deg3', 1.334, ...
    'C_deg4', 1.750, ...
    'N', 1.073, ...
    'O', 0.837);

fprintf('  alpha_units = %s\n', model.alpha_units);
fprintf('  thole_a     = %.3f\n', model.thole_a);

%% ------------------------------------------------------------------------
% 3. Build neutral periodic crystal system
% -------------------------------------------------------------------------

fprintf('\n[3] Building neutral crystal system...\n');

buildOpts = struct();
buildOpts.supercell_size = cfg.supercellSize;
buildOpts.bondScale = cfg.bondScale;
buildOpts.verbose = cfg.verbose;

tBuild = tic;
sys = builder.make_crystal_system(crystal, model, buildOpts);
buildTime = toc(tBuild);

% Ensure neutral/no assigned source charges.
if isfield(sys, 'site_charge')
    sys.site_charge(:) = 0;
end

if isfield(sys, 'site_is_active')
    sys.site_is_active(:) = false;
end

fprintf('\nBuilt neutral system:\n');
fprintf('  build time       = %.6f s\n', buildTime);
fprintf('  supercell_size   = [%d %d %d]\n', sys.supercell_size);
fprintf('  n_sites          = %d\n', sys.n_sites);
fprintf('  n_molecules      = %d\n', numel(sys.molecule_table.molecule_id));
fprintf('  n_complete       = %d\n', nnz(sys.molecule_table.is_complete_in_display));
fprintf('  total charge     = %+ .12e e\n', sum(sys.site_charge));
fprintf('  length unit      = %s\n', sys.units.length);
fprintf('  alpha unit       = %s\n', sys.units.alpha);

io.assert_atomic_units(sys);

lat0 = geom.get_lattice(sys);
Lmin0 = geom.shortest_lattice_translation(lat0.H);

fprintf('  lattice volume   = %.8e bohr^3\n', lat0.volume);
fprintf('  ||H*G - 2piI||   = %.3e\n', norm(lat0.H * lat0.G - 2*pi*eye(3), 'fro'));
fprintf('  Lmin             = %.8f bohr\n', Lmin0);
fprintf('  Lmin/2           = %.8f bohr\n', 0.5 * Lmin0);

if ~(cfg.ewald.rcut < 0.5 * Lmin0)
    error(['Periodic real-space cache requires rcut < Lmin/2.\n' ...
           '  rcut   = %.8f bohr\n' ...
           '  Lmin/2 = %.8f bohr\n' ...
           'Increase the supercell or reduce cfg.ewald.rcut.'], ...
           cfg.ewald.rcut, 0.5 * Lmin0);
end

%% ------------------------------------------------------------------------
% 4. Extract periodic polarization system
% -------------------------------------------------------------------------

fprintf('\n[4] Extracting periodic polarization system...\n');

polsys = builder.extract_polarization_system(sys, struct('mode', 'periodic'));
io.assert_atomic_units(polsys);

lat = geom.get_lattice(polsys);
V = lat.volume;

targetMask = logical(polsys.site_is_polarizable(:));
nPolSites = nnz(targetMask);

fprintf('  polsys.n_sites     = %d\n', polsys.n_sites);
fprintf('  polsys.is_periodic = %d\n', polsys.is_periodic);
fprintf('  polarizable sites  = %d\n', nPolSites);
fprintf('  total charge       = %+ .12e e\n', sum(polsys.site_charge));
fprintf('  lattice volume     = %.12e bohr^3\n', V);

if abs(sum(polsys.site_charge)) > 1e-12
    error('Dielectric tensor workflow expects a neutral system with zero fixed charges.');
end

%% ------------------------------------------------------------------------
% 5. Build zero-field problem and periodic P3M operator
% -------------------------------------------------------------------------

fprintf('\n[5] Building periodic P3M polarization operator...\n');

Ezero = zeros(polsys.n_sites, 3);

problem0 = thole.prepare_scf_problem(polsys, Ezero, cfg.scf);

tOp = tic;
op = local_build_periodic_p3m_operator(polsys, problem0, cfg);
opTime = toc(tOp);

fprintf('  operator build wall time = %.6f s\n', opTime);
fprintf('  operator mode            = %s\n', op.mode);
fprintf('  operator kind            = %s\n', op.kind);
fprintf('  operator backend         = %s\n', op.backend);
fprintf('  operator size            = %d x %d\n', op.size(1), op.size(2));

if isfield(op, 'info')
    if isfield(op.info, 'cache_time')
        fprintf('  cache time               = %.6f s\n', op.info.cache_time);
    end
    if isfield(op.info, 'real_cache_time')
        fprintf('  real cache time          = %.6f s\n', op.info.real_cache_time);
    end
    if isfield(op.info, 'p3m_cache_time')
        fprintf('  p3m cache time           = %.6f s\n', op.info.p3m_cache_time);
    end
    if isfield(op.info, 'nRealEntriesDirected')
        fprintf('  real directed entries    = %d\n', op.info.nRealEntriesDirected);
    end
    if isfield(op.info, 'nK')
        fprintf('  reciprocal k-vectors     = %d\n', op.info.nK);
    end
end

%% ------------------------------------------------------------------------
% 6. Solve response to uniform fields
% -------------------------------------------------------------------------

fprintf('\n[6] Solving response to uniform Cartesian fields...\n');

axesLabels = {'x', 'y', 'z'};
E0 = cfg.field_amplitude_au;

M_plus = zeros(3, 3);
M_minus = zeros(3, 3);
M_fd = zeros(3, 3);

P_plus = zeros(3, 3);
P_fd = zeros(3, 3);

solveTable = table();

for j = 1:3
    fprintf('\n  Field along %s: +%.6e a.u.\n', axesLabels{j}, E0);

    Efield = zeros(polsys.n_sites, 3);
    Efield(targetMask, j) = E0;

    problemPlus = thole.prepare_scf_problem(polsys, Efield, cfg.scf);

    tSolve = tic;
    [muPlus, infoPlus] = local_solve_selected(problemPlus, op, cfg);
    solveTimePlus = toc(tSolve);

    M_plus(:, j) = sum(muPlus, 1).';
    P_plus(:, j) = M_plus(:, j) ./ V;

    fprintf('    solve time       = %.6f s\n', solveTimePlus);
    fprintf('    converged        = %d\n', infoPlus.converged);
    fprintf('    relres           = %.12e\n', infoPlus.relres);
    fprintf('    ||mu||_F         = %.12e\n', norm(muPlus, 'fro'));
    fprintf('    M_ind / au       = [%+.8e %+.8e %+.8e]\n', ...
        M_plus(1,j), M_plus(2,j), M_plus(3,j));

    rowPlus = table();
    rowPlus.axis = string(axesLabels{j});
    rowPlus.sign = "+";
    rowPlus.solve_time_s = solveTimePlus;
    rowPlus.converged = logical(infoPlus.converged);
    rowPlus.relres = infoPlus.relres;
    rowPlus.mu_norm = norm(muPlus, 'fro');
    rowPlus.Mx = M_plus(1, j);
    rowPlus.My = M_plus(2, j);
    rowPlus.Mz = M_plus(3, j);

    solveTable = [solveTable; rowPlus]; %#ok<AGROW>

    if cfg.run_plus_minus_check
        fprintf('  Field along %s: -%.6e a.u.\n', axesLabels{j}, E0);

        EfieldMinus = zeros(polsys.n_sites, 3);
        EfieldMinus(targetMask, j) = -E0;

        problemMinus = thole.prepare_scf_problem(polsys, EfieldMinus, cfg.scf);

        tSolve = tic;
        [muMinus, infoMinus] = local_solve_selected(problemMinus, op, cfg);
        solveTimeMinus = toc(tSolve);

        M_minus(:, j) = sum(muMinus, 1).';
        M_fd(:, j) = 0.5 .* (M_plus(:, j) - M_minus(:, j));
        P_fd(:, j) = M_fd(:, j) ./ V;

        oddness = norm(M_plus(:,j) + M_minus(:,j)) / ...
            max(norm(M_plus(:,j) - M_minus(:,j)), eps);

        fprintf('    solve time       = %.6f s\n', solveTimeMinus);
        fprintf('    converged        = %d\n', infoMinus.converged);
        fprintf('    relres           = %.12e\n', infoMinus.relres);
        fprintf('    oddness check    = %.12e\n', oddness);
        fprintf('    M_ind / au       = [%+.8e %+.8e %+.8e]\n', ...
            M_minus(1,j), M_minus(2,j), M_minus(3,j));

        rowMinus = table();
        rowMinus.axis = string(axesLabels{j});
        rowMinus.sign = "-";
        rowMinus.solve_time_s = solveTimeMinus;
        rowMinus.converged = logical(infoMinus.converged);
        rowMinus.relres = infoMinus.relres;
        rowMinus.mu_norm = norm(muMinus, 'fro');
        rowMinus.Mx = M_minus(1, j);
        rowMinus.My = M_minus(2, j);
        rowMinus.Mz = M_minus(3, j);

        solveTable = [solveTable; rowMinus]; %#ok<AGROW>
    end
end

%% ------------------------------------------------------------------------
% 7. Compute chi and epsilon
% -------------------------------------------------------------------------

fprintf('\n[7] Computing susceptibility and dielectric tensor...\n');

if cfg.run_plus_minus_check
    M_response = M_fd;
    P_response = P_fd;
    responseLabel = 'central finite difference using +/- uniform fields';
else
    M_response = M_plus;
    P_response = P_plus;
    responseLabel = 'forward response using + uniform fields';
end

chi = P_response ./ E0;
epsilon = eye(3) + 4*pi*chi;

% -------------------------------------------------------------------------
% Alternative dielectric conventions / local-field diagnostics
% -------------------------------------------------------------------------

I3 = eye(3);

% 1. Direct response already computed:
epsilon_direct = epsilon;

% 2. Isotropic Lorentz local-field correction applied to the direct
% susceptibility.
%
% Model:
%   P = chi_direct * E_local
%   E_local = E_macro + L * P
%   L = 4*pi/3 I
%
% Therefore:
%   P = chi_direct * (E_macro + L P)
%   (I - chi_direct L) P = chi_direct E_macro
%   chi_macro = (I - chi_direct L)^(-1) chi_direct
%
% Warning:
%   This may over-correct if chi_direct already includes part of the
%   microscopic local-field response. Treat as diagnostic.
L_lorentz = (4*pi/3) * I3;
chi_lorentz_from_direct = (I3 - chi * L_lorentz) \ chi;
epsilon_lorentz_from_direct = I3 + 4*pi * chi_lorentz_from_direct;

epsilon_lorentz_sym = 0.5 * (epsilon_lorentz_from_direct + epsilon_lorentz_from_direct.');

[vecLorentz, valLorentz] = eig(epsilon_lorentz_sym);
epsLorentzPrincipal = diag(valLorentz);
[epsLorentzPrincipal, order] = sort(epsLorentzPrincipal, 'ascend');
vecLorentz = vecLorentz(:, order);

fprintf('\nAlternative dielectric conventions:\n');

fprintf('\nDirect periodic-field epsilon:\n');
disp(epsilon_direct);

fprintf('Direct periodic-field principal epsilon:\n');
disp(sort(eig(0.5*(epsilon_direct + epsilon_direct.')), 'ascend').');

fprintf('\nLorentz-corrected-from-direct epsilon:\n');
disp(epsilon_lorentz_from_direct);

fprintf('Lorentz-corrected-from-direct principal epsilon:\n');
disp(epsLorentzPrincipal.');

fprintf('Lorentz-corrected-from-direct principal n:\n');
disp(sqrt(epsLorentzPrincipal).');

epsilonSym = 0.5 .* (epsilon + epsilon.');
epsilonAnti = 0.5 .* (epsilon - epsilon.');

[eVec, eValMat] = eig(epsilonSym);
epsPrincipal = diag(eValMat);

fprintf('  response convention = %s\n', responseLabel);
fprintf('  volume              = %.12e bohr^3\n', V);
fprintf('  field amplitude     = %.12e a.u.\n', E0);

fprintf('\nInduced cell dipole response M columns / au:\n');
disp(M_response);

fprintf('Susceptibility tensor chi, Cartesian:\n');
disp(chi);

fprintf('Dielectric tensor epsilon = I + 4*pi*chi, Cartesian:\n');
disp(epsilon);

fprintf('Symmetrized dielectric tensor:\n');
disp(epsilonSym);

fprintf('Antisymmetric part norm / Frobenius:\n');
fprintf('  ||anti(eps)||_F = %.12e\n', norm(epsilonAnti, 'fro'));
fprintf('  ||sym(eps)||_F  = %.12e\n', norm(epsilonSym, 'fro'));
fprintf('  ratio           = %.12e\n', ...
    norm(epsilonAnti, 'fro') / max(norm(epsilonSym, 'fro'), eps));

fprintf('\nPrincipal dielectric values from symmetrized tensor:\n');
fprintf('  eps_principal = [%.8f %.8f %.8f]\n', ...
    epsPrincipal(1), epsPrincipal(2), epsPrincipal(3));

fprintf('\nPrincipal axes columns, Cartesian:\n');
disp(eVec);

%% ------------------------------------------------------------------------
% 8. Summary table
% -------------------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('Dielectric tensor workflow summary\n');
fprintf('============================================================\n');

disp(solveTable);

summary = table();
summary.supercell = string(sprintf('[%d %d %d]', cfg.supercellSize));
summary.nSites = polsys.n_sites;
summary.nPolSites = nPolSites;
summary.volume_bohr3 = V;
summary.field_amplitude_au = E0;
summary.boundary = string(cfg.ewald.boundary);
summary.mesh = string(sprintf('[%d %d %d]', cfg.p3m.mesh_size));
summary.eps_xx = epsilon(1,1);
summary.eps_yy = epsilon(2,2);
summary.eps_zz = epsilon(3,3);
summary.eps_xy = epsilon(1,2);
summary.eps_xz = epsilon(1,3);
summary.eps_yz = epsilon(2,3);
summary.eps_principal_1 = epsPrincipal(1);
summary.eps_principal_2 = epsPrincipal(2);
summary.eps_principal_3 = epsPrincipal(3);
summary.antisym_norm = norm(epsilonAnti, 'fro');

fprintf('\nCompact summary:\n');
disp(summary);

fprintf('\nWorkflow completed successfully.\n');

%% =========================================================================
% Local helpers
% =========================================================================

function op = local_build_periodic_p3m_operator(polsys, problem, cfg)
op = thole.make_polarization_operator(polsys, problem, ...
    'Mode', 'periodic_p3m', ...
    'Solver', cfg.solver.method, ...
    'Backend', cfg.operator.backend, ...
    'UseThole', cfg.operator.use_thole, ...
    'Softening', cfg.operator.softening, ...
    'Rcut', cfg.ewald.rcut, ...
    'Alpha', cfg.ewald.alpha, ...
    'Kcut', cfg.ewald.kcut, ...
    'Boundary', cfg.ewald.boundary, ...
    'MeshSize', cfg.p3m.mesh_size, ...
    'AssignmentOrder', cfg.p3m.assignment_order, ...
    'DerivativeMode', cfg.p3m.derivative_mode, ...
    'InfluenceMode', cfg.p3m.influence_mode, ...
    'FDStencil', cfg.p3m.fd_stencil, ...
    'DeconvolveAssignment', cfg.p3m.deconvolve_assignment, ...
    'DeconvolutionFloor', cfg.p3m.deconvolution_floor, ...
    'AliasRange', cfg.p3m.alias_range, ...
    'UseMex', cfg.operator.use_mex, ...
    'Profile', cfg.operator.profile, ...
    'Verbose', cfg.operator.verbose);
end

function [mu, solveInfo] = local_solve_selected(problem, op, cfg)
switch lower(cfg.solver.method)
    case 'direct'
        solveOpts = struct();
        solveOpts.compute_residual = cfg.solver.compute_residual;

        [mu, solveInfo] = thole.solve_scf_direct(problem, op, solveOpts);

    case 'jacobi'
        solveOpts = struct();
        solveOpts.tol = cfg.scf.tol;
        solveOpts.max_iter = cfg.scf.maxIter;
        solveOpts.mixing = 0.6;
        solveOpts.verbose = cfg.scf.verbose;

        [mu, solveInfo] = thole.solve_scf_jacobi(problem, op, solveOpts);

    case 'gmres'
        solveOpts = struct();
        solveOpts.tol = cfg.scf.tol;
        solveOpts.max_iter = cfg.scf.maxIter;
        solveOpts.restart = cfg.solver.gmres_restart;
        solveOpts.verbose = cfg.scf.verbose;

        [mu, solveInfo] = thole.solve_scf_gmres(problem, op, solveOpts);

    case 'sor'
        solveOpts = struct();
        solveOpts.tol = cfg.scf.tol;
        solveOpts.max_iter = cfg.scf.maxIter;
        solveOpts.omega = 0.52;
        solveOpts.stop_metric = 'max_dmu';
        solveOpts.residual_every = 25;
        solveOpts.verbose = cfg.scf.verbose;

        [mu, solveInfo] = thole.solve_scf_sor(problem, op, solveOpts);

    otherwise
        error('Unsupported cfg.solver.method "%s".', cfg.solver.method);
end
end