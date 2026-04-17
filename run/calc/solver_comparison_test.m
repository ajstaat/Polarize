%% run_test_nonperiodic_solver_comparison
clear; clc; close all;

fprintf('=== nonperiodic solver comparison test ===\n');

rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
filename   = fullfile(rootFolder, 'a_0.0_CONTCAR.vasp');

%% ------------------------------------------------------------------------
% User controls
supercellSize = [2 5 1];
bondScale = 1.20;

% Optional nonperiodic cutoff in bohr.
% Set [] for exact all-pairs.
rcut_bohr = 15.0; %[];

% Charged-pair assignment
pairCharges = [+1 -1];

% Solver controls
sorOmega = 0.90;
jacobiMixing = 0.60;      % plain Jacobi
jacobiMaxIter = 200;

%% ------------------------------------------------------------------------
% Import crystal template
crystal = io.import_contcar_as_crystal(filename, ...
    'BondScale', 1.20, ...
    'WrapFractional', true, ...
    'SortMolecules', false);

fprintf('Imported crystal template:\n');
fprintf('  nSites     = %d\n', size(crystal.cart_coords, 1));
fprintf('  nBaseMols  = %d\n', numel(unique(crystal.base_mol_id)));

%% ------------------------------------------------------------------------
% Model with updated class-based polarizabilities
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
% Build working supercell
opts = struct();
opts.supercell_size = supercellSize;
opts.bondScale = bondScale;
opts.verbose = true;
opts.removeMolIDs = [];
opts.activeMolIDs = [];

sys = builder.make_crystal_system(crystal, model, opts);

fprintf('\nBuilt working system:\n');
fprintf('  n_unit_sites = %d\n', sys.n_unit_sites);
fprintf('  n_cells      = %d\n', sys.n_cells);
fprintf('  n_sites      = %d\n', sys.n_sites);

%% ------------------------------------------------------------------------
% Automatically choose a complete center reference molecule and its
% same-stack first-shell neighbor
[refID, refSummary] = builder.choose_center_reference_molecule(sys, ...
    'RequireComplete', true, ...
    'Verbose', true);

desc = builder.complete_molecule_descriptors_relative_to_reference(sys, refID, ...
    'Verbose', true);

sameStack = builder.select_same_stack_neighbor(desc, 1, ...
    'Direction', 'either', ...
    'Verbose', true);

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

if height(sameStack.match_table) > 1
    fprintf('  same-stack shell-1 candidates = %s\n', ...
        mat2str(sameStack.match_table.molecule_id(:).'));
end

%% ------------------------------------------------------------------------
% Apply uniform cation/anion charges and disable polarizability on those sites
sys = builder.apply_molecule_charges(sys, [refID nbrID], ...
    'Mode', 'uniform', ...
    'TotalCharges', pairCharges, ...
    'SetActive', true, ...
    'DisablePolarizabilityOnCharged', true, ...
    'ZeroExistingCharges', true, ...
    'Verbose', true);

%% ------------------------------------------------------------------------
% Basic charged-pair diagnostics
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
fprintf('  active sites total              = %d\n', nnz(sys.site_is_active));

if abs(sum(sys.site_charge(refIdx)) - pairCharges(1)) > 1e-12
    error('Reference molecule charge assignment failed.');
end

if abs(sum(sys.site_charge(nbrIdx)) - pairCharges(2)) > 1e-12
    error('Neighbor molecule charge assignment failed.');
end

if any(sys.site_is_polarizable(refIdx))
    error('Reference molecule still has polarizable sites enabled.');
end

if any(sys.site_is_polarizable(nbrIdx))
    error('Neighbor molecule still has polarizable sites enabled.');
end

if ~all(sys.site_is_active(refIdx)) || ~all(sys.site_is_active(nbrIdx))
    error('Active mask was not set correctly on charged molecules.');
end

%% ------------------------------------------------------------------------
% Canonical AU system
params = struct();
params.ewald = struct();
params.ewald.mode = 'nonperiodic';

polsys = builder.extract_polarization_system(sys, params);
io.assert_atomic_units(polsys);

fprintf('\nCanonical polarization system:\n');
fprintf('  n_sites             = %d\n', polsys.n_sites);
fprintf('  n_polarizable_sites = %d\n', nnz(polsys.site_is_polarizable));

%% ------------------------------------------------------------------------
% External field
fieldParams = struct();
fieldParams.field = struct();
fieldParams.field.include_external_charges = true;
fieldParams.field.exclude_self = true;
fieldParams.field.softening = 0.0;

if ~isempty(rcut_bohr)
    fieldParams.field.rcut = rcut_bohr;
end

Eext = calc.compute_external_field(polsys, fieldParams);

if isempty(Eext) || ~isequal(size(Eext), [polsys.n_sites, 3])
    error('External field was not built correctly.');
end

if any(~isfinite(Eext(:)))
    error('External field contains non-finite values.');
end

%% ------------------------------------------------------------------------
% Quiet settings for operator assembly
scfAssembly = struct();
scfAssembly.use_thole = true;
scfAssembly.softening = 0.0;
scfAssembly.tol = 1e-8;
scfAssembly.maxIter = 500;
scfAssembly.verbose = true;        % suppress Tpol row-by-row printing
scfAssembly.printEvery = 250;
scfAssembly.residualEvery = 250;
scfAssembly.stopMetric = 'relres';

if ~isempty(rcut_bohr)
    scfAssembly.rcut = rcut_bohr;
end

problemAssembly = thole.prepare_scf_problem(polsys, Eext, scfAssembly);
[Tpol, opinfo] = ewald.assemble_nonperiodic_interaction_matrix(polsys, problemAssembly, struct(), scfAssembly);

fprintf('\nOperator diagnostics:\n');
fprintf('  active-space operator size = %d x %d\n', size(Tpol,1), size(Tpol,2));
fprintf('  nPolSites                  = %d\n', problemAssembly.nPolSites);
fprintf('  opinfo.space               = %s\n', opinfo.space);

if isfield(opinfo, 'rcut') && isfinite(opinfo.rcut)
    nPol = problemAssembly.nPolSites;
    nPairsFull = nPol * (nPol - 1) / 2;
    fprintf('  opinfo.rcut                = %.6f bohr\n', opinfo.rcut);
    fprintf('  pair blocks kept           = %d / %d (%.4f%%)\n', ...
        opinfo.nPairBlocksKept, nPairsFull, ...
        100 * opinfo.nPairBlocksKept / nPairsFull);
else
    fprintf('  opinfo.rcut                = exact (no cutoff)\n');
end

if ~isequal(size(Tpol), [3*problemAssembly.nPolSites, 3*problemAssembly.nPolSites])
    error('Tpol size mismatch.');
end

if any(~isfinite(Tpol(:)))
    error('Tpol contains non-finite values.');
end

%% ------------------------------------------------------------------------
% Direct solve
fprintf('\n--- direct solve ---\n');
problemDirect = thole.prepare_scf_problem(polsys, Eext, scfAssembly);

tic;
[mu_direct, direct] = thole.solve_scf_direct(problemDirect, Tpol);
time_direct = toc;

E_direct = calc.compute_total_energy_active_space(polsys, problemDirect, mu_direct, Eext, Tpol);
E_direct_eV = E_direct.total * 27.211386245988;
relres_direct = thole.compute_active_space_relres(problemDirect, Tpol, mu_direct);

fprintf('direct: relres = %.3e | E = %+0.12f Ha | %+0.12f eV | time = %.3f s\n', ...
    relres_direct, E_direct.total, E_direct_eV, time_direct);

%% ------------------------------------------------------------------------
% SOR solve
fprintf('\n--- SOR solve ---\n');
scfSOR = scfAssembly;
scfSOR.verbose = true;
scfSOR.printEvery = 5;
scfSOR.residualEvery = 5;
scfSOR.omega = sorOmega;

problemSOR = thole.prepare_scf_problem(polsys, Eext, scfSOR);

tic;
[mu_sor, scf_sor] = thole.solve_scf_sor(problemSOR, Tpol);
time_sor = toc;

E_sor = calc.compute_total_energy_active_space(polsys, problemSOR, mu_sor, Eext, Tpol);
E_sor_eV = E_sor.total * 27.211386245988;
relres_sor = thole.compute_active_space_relres(problemSOR, Tpol, mu_sor);

fprintf('sor:    relres = %.3e | E = %+0.12f Ha | %+0.12f eV | iters = %d | time = %.3f s\n', ...
    relres_sor, E_sor.total, E_sor_eV, scf_sor.nIter, time_sor);

%% ------------------------------------------------------------------------
% Jacobi solve
% fprintf('\n--- Jacobi solve ---\n');
% scfJac = scfAssembly;
% scfJac.verbose = true;
% scfJac.printEvery = 5;       % print every Jacobi iteration
% scfJac.residualEvery = 5;    % compute relres every Jacobi iteration
% scfJac.mixing = jacobiMixing;
% scfJac.maxIter = jacobiMaxIter;
% 
% problemJac = thole.prepare_scf_problem(polsys, Eext, scfJac);
% 
% tic;
% [mu_jac, scf_jac] = thole.solve_scf_iterative(problemJac, Tpol);
% time_jac = toc;
% 
% E_jac = calc.compute_total_energy_active_space(polsys, problemJac, mu_jac, Eext, Tpol);
% E_jac_eV = E_jac.total * 27.211386245988;
% relres_jac = thole.compute_active_space_relres(problemJac, Tpol, mu_jac);
% 
% fprintf('jacobi: relres = %.3e | E = %+0.12f Ha | %+0.12f eV | iters = %d | time = %.3f s | converged = %d\n', ...
%     relres_jac, E_jac.total, E_jac_eV, scf_jac.nIter, time_jac, scf_jac.converged);

%% ------------------------------------------------------------------------
% Matrix-free iterative solve
fprintf('\n--- matrix-free iterative solve ---\n');
scfJac = scfAssembly;
scfJac.verbose = true;
scfJac.printEvery = 1;
scfJac.residualEvery = 1;
scfJac.mixing = jacobiMixing;
scfJac.maxIter = jacobiMaxIter;

tic;
[mu_jac, scf_jac] = thole.solve_scf_iterative(polsys, Eext, scfJac);
time_jac = toc;

% Still build problem for apples-to-apples reporting with explicit Tpol
problemJac = thole.prepare_scf_problem(polsys, Eext, scfJac);
E_jac = calc.compute_total_energy_active_space(polsys, problemJac, mu_jac, Eext, Tpol);
E_jac_eV = E_jac.total * 27.211386245988;
relres_jac = thole.compute_active_space_relres(problemJac, Tpol, mu_jac);

fprintf('iter: relres = %.3e | E = %+0.12f Ha | %+0.12f eV | iters = %d | time = %.3f s | converged = %d\n', ...
    relres_jac, E_jac.total, E_jac_eV, scf_jac.nIter, time_jac, scf_jac.converged);

%% ------------------------------------------------------------------------
% Compare solutions
fprintf('\n--- solver comparisons ---\n');

mu_diff_direct_sor = max(sqrt(sum((mu_direct - mu_sor).^2, 2)));
mu_diff_direct_jac = max(sqrt(sum((mu_direct - mu_jac).^2, 2)));
mu_diff_sor_jac    = max(sqrt(sum((mu_sor - mu_jac).^2, 2)));

fprintf('max |mu_direct - mu_sor| = %.16e\n', mu_diff_direct_sor);
fprintf('max |mu_direct - mu_jac| = %.16e\n', mu_diff_direct_jac);
fprintf('max |mu_sor    - mu_jac| = %.16e\n', mu_diff_sor_jac);

fprintf('|E_direct - E_sor| = %.16e Ha  (%.16e eV)\n', ...
    abs(E_direct.total - E_sor.total), abs(E_direct_eV - E_sor_eV));
fprintf('|E_direct - E_jac| = %.16e Ha  (%.16e eV)\n', ...
    abs(E_direct.total - E_jac.total), abs(E_direct_eV - E_jac_eV));
fprintf('|E_sor    - E_jac| = %.16e Ha  (%.16e eV)\n', ...
    abs(E_sor.total - E_jac.total), abs(E_sor_eV - E_jac_eV));

%% ------------------------------------------------------------------------
% Assertions
if ~isfinite(direct.relres) || direct.relres > 1e-10
    error('Direct solver relative residual too large: %.3e', direct.relres);
end

if ~scf_sor.converged
    error('SOR solver did not converge.');
end

if ~isfinite(relres_sor) || relres_sor > 1e-8
    error('SOR relative residual too large: %.3e', relres_sor);
end

if mu_diff_direct_sor > 1e-8
    error('Direct and SOR dipoles disagree: %.3e', mu_diff_direct_sor);
end

if abs(E_direct.total - E_sor.total) > 1e-8
    error('Direct and SOR energies disagree: %.3e Ha', abs(E_direct.total - E_sor.total));
end

% Jacobi is diagnostic here, not a hard-pass requirement.
if ~scf_jac.converged
    fprintf('\nJacobi did not converge for this system.\n');
end

%% ------------------------------------------------------------------------
% Optional plot: SOR dipoles
BOHR2ANG = 0.529177210903;

fig2 = figure(2);
clf(fig2);
fig2.Color = 'w';
ax2 = axes('Parent', fig2);
hold(ax2, 'on');
grid(ax2, 'on');
box(ax2, 'on');
view(ax2, 3);
axis(ax2, 'equal');
xlabel(ax2, 'x (Å)');
ylabel(ax2, 'y (Å)');
zlabel(ax2, 'z (Å)');

Xplot = sys.site_pos;
polIdx = find(sys.site_is_polarizable);

plotSkip = 1;
sel = 1:plotSkip:numel(polIdx);
plotIdx = polIdx(sel);

arrowScale = 50.0;

plot3(ax2, Xplot(:,1), Xplot(:,2), Xplot(:,3), '.', ...
    'Color', [0.75 0.75 0.75], 'MarkerSize', 4);

scatter3(ax2, Xplot(refIdx,1), Xplot(refIdx,2), Xplot(refIdx,3), 30, 'filled');
scatter3(ax2, Xplot(nbrIdx,1), Xplot(nbrIdx,2), Xplot(nbrIdx,3), 30, 'filled');

quiver3(ax2, ...
    Xplot(plotIdx,1), Xplot(plotIdx,2), Xplot(plotIdx,3), ...
    arrowScale * mu_sor(plotIdx,1) * BOHR2ANG, ...
    arrowScale * mu_sor(plotIdx,2) * BOHR2ANG, ...
    arrowScale * mu_sor(plotIdx,3) * BOHR2ANG, ...
    0, 'k', 'LineWidth', 1.0, 'MaxHeadSize', 0.6);

title(ax2, sprintf(['Converged induced dipoles (SOR) | ref %d (+1), nbr %d (-1) | ' ...
    'supercell [%d %d %d]'], ...
    refID, nbrID, supercellSize));

%% ------------------------------------------------------------------------
% Final report
fprintf('\nDONE: nonperiodic solver comparison completed.\n');
fprintf('  supercell size   = [%d %d %d]\n', supercellSize);
fprintf('  ref ID           = %d\n', refID);
fprintf('  nbr ID           = %d\n', nbrID);
fprintf('  direct energy    = %+0.12f Ha  | %+0.12f eV\n', E_direct.total, E_direct_eV);
fprintf('  sor energy       = %+0.12f Ha  | %+0.12f eV\n', E_sor.total, E_sor_eV);
fprintf('  jacobi energy    = %+0.12f Ha  | %+0.12f eV\n', E_jac.total, E_jac_eV);
fprintf('  direct relres    = %.3e\n', relres_direct);
fprintf('  sor relres       = %.3e\n', relres_sor);
fprintf('  jacobi relres    = %.3e\n', relres_jac);
fprintf('  sor nIter        = %d\n', scf_sor.nIter);
fprintf('  jacobi nIter     = %d\n', scf_jac.nIter);
if ~isempty(rcut_bohr)
    fprintf('  cutoff           = %.6f bohr\n', rcut_bohr);
else
    fprintf('  cutoff           = exact (no cutoff)\n');
end