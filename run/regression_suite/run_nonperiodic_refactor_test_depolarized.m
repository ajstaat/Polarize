function run_nonperiodic_refactor_test_depolarized
clc;

fprintf('=== nonperiodic refactor regression test: depolarized active molecules ===\n');

[crystal, model, opts0, ~, molA, molB] = make_standard_nonperiodic_case([2 2 1]);

opts = opts0;
opts.activeMolIDs = [molA, molB];
opts.depolarizeActiveMolecules = true;

params = util.default_params();
params.scf.use_thole = true;
params.scf.tol = 1e-12;
params.scf.maxIter = 1000;
params.scf.mixing = 0.35;
params.scf.omega = 1.0;
params.output.computeEnergy = true;

params_iter = params; params_iter.scf.solver = 'matrix_iterative';
params_dir  = params; params_dir.scf.solver  = 'direct';
params_sor  = params; params_sor.scf.solver  = 'sor';

res_iter = calc.run_polarization_calc(crystal, model, opts, params_iter);
res_dir  = calc.run_polarization_calc(crystal, model, opts, params_dir);
res_sor  = calc.run_polarization_calc(crystal, model, opts, params_sor);

err_iter = max(sqrt(sum((res_iter.mu - res_dir.mu).^2, 2)));
err_sor  = max(sqrt(sum((res_sor.mu  - res_dir.mu).^2, 2)));

fprintf('max |mu_iter - mu_direct| = %.3e\n', err_iter);
fprintf('max |mu_sor  - mu_direct| = %.3e\n', err_sor);

assert(err_iter < 1e-8);
assert(err_sor  < 1e-8);

fprintf('PASS\n');
end