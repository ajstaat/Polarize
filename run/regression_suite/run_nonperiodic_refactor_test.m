function run_nonperiodic_refactor_test
clc;

fprintf('=== nonperiodic refactor regression test ===\n');

[crystal, model, opts0, ~, molA, molB] = make_standard_nonperiodic_case([2 2 1]);

opts = opts0;
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

[sys, polsys, Eext, problem, Tpol] = build_nonperiodic_test_problem(crystal, model, opts, params);

fprintf('nSites = %d, nPolSites = %d\n', problem.nSites, problem.nPolSites);

fprintf('\n--- operator mapping ---\n');
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

fprintf('\n--- direct solve ---\n');
[mu_dir, ~] = thole.solve_scf_direct(problem, Tpol);
fprintf('direct solve complete\n');

fprintf('\n--- matrix_iterative solve ---\n');
[mu_mat, scf_mat] = thole.solve_scf_matrix_iterative(problem, Tpol);
err_mat_dir = max(sqrt(sum((mu_mat - mu_dir).^2, 2)));
fprintf('matrix_iterative converged=%d iter=%d\n', scf_mat.converged, scf_mat.nIter);
fprintf('max |mu_matrix - mu_direct| = %.16e\n', err_mat_dir);

fprintf('\n--- sor solve ---\n');
problem_sor = thole.prepare_scf_problem(polsys, Eext, params.scf);
[mu_sor, scf_sor] = thole.solve_scf_sor(problem_sor, Tpol);
err_sor_dir = max(sqrt(sum((mu_sor - mu_dir).^2, 2)));
fprintf('sor converged=%d iter=%d\n', scf_sor.converged, scf_sor.nIter);
fprintf('max |mu_sor - mu_direct| = %.16e\n', err_sor_dir);

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