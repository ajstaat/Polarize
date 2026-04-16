function run_nonperiodic_refactor_test_initial_mu
clc;

fprintf('=== nonperiodic refactor regression test: nonzero initial guess ===\n');

[crystal, model, opts0, params0, molA, molB] = make_standard_nonperiodic_case([2 2 1]);

opts = opts0;
opts.activeMolIDs = [molA, molB];

sys = builder.make_crystal_system(crystal, model, opts);
sys = builder.assign_point_charges(sys, model.charge_pattern);

rng(1);
initial_mu = 1e-3 * randn(sys.n_sites, 3);

params = util.default_params();
params.scf.use_thole = true;
params.scf.tol = 1e-12;
params.scf.maxIter = 1000;
params.scf.mixing = 0.35;
params.scf.omega = 1.0;
params.scf.initial_mu = initial_mu;
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