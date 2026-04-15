function run_nonperiodic_refactor_test_depolarized
clc;

fprintf('=== nonperiodic refactor regression test: depolarized active molecules ===\n');

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
model.polarizable_types = {'C','N'};
model.alpha_by_type = struct();
model.alpha_by_type.C = 10.0;
model.alpha_by_type.N = 8.0;
model.thole_a = 0.39;

model.charge_pattern = struct();
model.charge_pattern.site_label = {'C1','H1'};
model.charge_pattern.delta_q = [0.10, -0.10];

opts0 = struct();
opts0.supercell_size = [2 2 1];
opts0.verbose = false;
opts0.removeMolIDs = [];
opts0.activeMolIDs = [];

params0 = util.default_params();
params0.scf.solver = 'matrix_iterative';

result0 = calc.run_polarization_calc(crystal, model, opts0, params0);
sys0 = result0.sys;

molA = builder.find_unique_molecule_id(sys0, 1, [0 0 0]);
molB = builder.find_unique_molecule_id(sys0, 1, [1 0 0]);

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

params_iter = params;
params_iter.scf.solver = 'matrix_iterative';

params_dir = params;
params_dir.scf.solver = 'direct';

params_sor = params;
params_sor.scf.solver = 'sor';

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