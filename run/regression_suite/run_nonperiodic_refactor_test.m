function run_nonperiodic_refactor_test
clc;

fprintf('=== nonperiodic refactor regression test ===\n');

% -------------------------------------------------------------------------
% Small inline case patterned after your existing run scripts
% -------------------------------------------------------------------------
crystal = struct();
crystal.cellpar = [12.3, 14.1, 9.8, 90.0, 101.2, 88.4];
crystal.lattice = [];
crystal.frac_coords = [
    0.10 0.20 0.30
    0.15 0.22 0.35
    0.60 0.70 0.80
];
crystal.cart_coords = [];
crystal.mol_id = [1; 1; 2];
crystal.site_label = {'C1'; 'H1'; 'N1'};
crystal.site_type = {'C'; 'H'; 'N'};

model = struct();
model.polarizable_types = {'C', 'N'};
model.alpha_by_type = struct();
model.alpha_by_type.C = 10.0;
model.alpha_by_type.N = 8.0;
model.thole_a = 0.39;

model.charge_pattern = struct();
model.charge_pattern.site_label = {'C1', 'H1'};
model.charge_pattern.delta_q = [0.10, -0.10];

% -------------------------------------------------------------------------
% Neutral build for molecule lookup
% -------------------------------------------------------------------------
opts0 = struct();
opts0.supercell_size = [2 2 1];
opts0.verbose = false;
opts0.removeMolIDs = [];
opts0.activeMolIDs = [];

params0 = util.default_params();
params0.scf.solver = 'matrix_iterative';
params0.scf.tol = 1e-10;
params0.scf.maxIter = 500;
params0.scf.mixing = 0.35;
params0.scf.use_thole = true;

result0 = calc.run_polarization_calc(crystal, model, opts0, params0);
sys0 = result0.sys;

molA = builder.find_unique_molecule_id(sys0, 1, [0 0 0]);
molB = builder.find_unique_molecule_id(sys0, 1, [1 0 0]);

opts = struct();
opts.supercell_size = [2 2 1];
opts.verbose = false;
opts.removeMolIDs = [];
opts.activeMolIDs = [molA, molB];

params = util.default_params();
params.scf.use_thole = true;
params.scf.tol = 1e-12;
params.scf.maxIter = 1000;
params.scf.mixing = 0.35;
params.scf.omega = 1.0;
params.output.computeEnergy = true;
params.output.computeEnergyByMolecule = true;
params.output.computeDipoleDipoleDecomp = true;

% -------------------------------------------------------------------------
% Build sys + Eext once
% -------------------------------------------------------------------------
sys = builder.make_crystal_system(crystal, model, opts);
sys = builder.assign_point_charges(sys, model.charge_pattern);
Eext = calc.compute_external_field(sys, params);
problem = thole.prepare_scf_problem(sys, Eext, params.scf);

fprintf('nSites = %d, nPolSites = %d\n', problem.nSites, problem.nPolSites);

% -------------------------------------------------------------------------
% 1) Operator mapping test
% -------------------------------------------------------------------------
fprintf('\n--- operator mapping ---\n');

[Tpol, ~] = ewald.assemble_nonperiodic_interaction_matrix(sys, problem, params.ewald, params.scf);
Tfull = ewald.expand_active_operator_to_full(problem, Tpol);

idx = problem.activeVecIdx;
map_err = norm(Tfull(idx, idx) - Tpol, 'fro');
fprintf('||Tfull(active,active) - Tpol||_F = %.3e\n', map_err);
assert(map_err < 1e-14, 'Active-space expansion mismatch.');

mask = true(3*problem.nSites, 3*problem.nSites);
mask(idx, idx) = false;
outside_err = norm(Tfull(mask));
fprintf('||Tfull(outside active-active)|| = %.3e\n', outside_err);
assert(outside_err < 1e-14, 'Expanded T has nonzero entries outside active-active block.');

% -------------------------------------------------------------------------
% 2) Direct solve
% -------------------------------------------------------------------------
fprintf('\n--- direct solve ---\n');

[mu_dir, ~] = thole.solve_scf_direct(problem, Tpol);
fprintf('direct solve complete\n');

% -------------------------------------------------------------------------
% 3) Matrix-iterative solve
% -------------------------------------------------------------------------
fprintf('\n--- matrix_iterative solve ---\n');

[mu_mat, scf_mat] = thole.solve_scf_matrix_iterative(problem, Tpol);
err_mat_dir = max(sqrt(sum((mu_mat - mu_dir).^2, 2)));

fprintf('matrix_iterative converged=%d iter=%d\n', scf_mat.converged, scf_mat.nIter);
fprintf('max |mu_matrix - mu_direct| = %.16e\n', err_mat_dir);

% -------------------------------------------------------------------------
% 4) SOR solve
% -------------------------------------------------------------------------
fprintf('\n--- sor solve ---\n');

problem_sor = thole.prepare_scf_problem(sys, Eext, params.scf);
[mu_sor, scf_sor] = thole.solve_scf_sor(problem_sor, Tpol);
err_sor_dir = max(sqrt(sum((mu_sor - mu_dir).^2, 2)));

fprintf('sor converged=%d iter=%d\n', scf_sor.converged, scf_sor.nIter);
fprintf('max |mu_sor - mu_direct| = %.16e\n', err_sor_dir);

% -------------------------------------------------------------------------
% 5) Residual checks
% -------------------------------------------------------------------------
fprintf('\n--- residual checks ---\n');

relres_dir = thole.compute_active_space_relres(problem, Tpol, mu_dir);
relres_mat = thole.compute_active_space_relres(problem, Tpol, mu_mat);
relres_sor = thole.compute_active_space_relres(problem, Tpol, mu_sor);

fprintf('relres(direct) = %.3e\n', relres_dir);
fprintf('relres(matrix) = %.3e\n', relres_mat);
fprintf('relres(sor)    = %.3e\n', relres_sor);

assert(relres_dir < 1e-10, 'Direct residual too large.');
assert(relres_mat < 1e-8,  'Matrix-iterative residual too large.');
assert(relres_sor < 1e-8,  'SOR residual too large.');

% -------------------------------------------------------------------------
% 6) End-to-end driver comparisons
% -------------------------------------------------------------------------
fprintf('\n--- run_polarization_calc end-to-end ---\n');

params_iter = params;
params_iter.scf.solver = 'matrix_iterative';

params_dir = params;
params_dir.scf.solver = 'direct';

params_sor = params;
params_sor.scf.solver = 'sor';
params_sor.scf.omega = 1.0;

res_iter = calc.run_polarization_calc(crystal, model, opts, params_iter);
res_dir  = calc.run_polarization_calc(crystal, model, opts, params_dir);
res_sor  = calc.run_polarization_calc(crystal, model, opts, params_sor);

err_iter_dir = max(sqrt(sum((res_iter.mu - res_dir.mu).^2, 2)));
err_sor_dir  = max(sqrt(sum((res_sor.mu  - res_dir.mu ).^2, 2)));

fprintf('driver: max |mu_iter - mu_direct| = %.16e\n', err_iter_dir);
fprintf('driver: max |mu_sor  - mu_direct| = %.16e\n', err_sor_dir);

if isfield(res_iter, 'energy') && isfield(res_dir, 'energy') && ...
        isfield(res_iter.energy, 'total') && isfield(res_dir.energy, 'total')
    fprintf('E_total(iter) = %.16e\n', res_iter.energy.total);
    fprintf('E_total(dir ) = %.16e\n', res_dir.energy.total);
    fprintf('E_total(sor ) = %.16e\n', res_sor.energy.total);
    fprintf('|E_iter - E_dir| = %.16e\n', abs(res_iter.energy.total - res_dir.energy.total));
    fprintf('|E_sor  - E_dir| = %.16e\n', abs(res_sor.energy.total - res_dir.energy.total));
end

assert(err_iter_dir < 1e-8, 'Driver matrix_iterative disagrees with direct.');
assert(err_sor_dir  < 1e-8, 'Driver SOR disagrees with direct.');

fprintf('\nPASS: nonperiodic refactor regression test completed.\n');
end