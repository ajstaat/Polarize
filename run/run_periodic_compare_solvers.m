clear; clc;

% -----------------------------
% Define crystal
% -----------------------------
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

% -----------------------------
% Define model
% -----------------------------
model = struct();
model.polarizable_types = {'C', 'N'};
model.alpha_by_type = struct();
model.alpha_by_type.C = 10.0;
model.alpha_by_type.N = 8.0;
model.thole_a = 0.39;

model.charge_pattern = struct();
model.charge_pattern.site_label = {'C1', 'H1'};
model.charge_pattern.delta_q = [0.10, -0.10];

% -----------------------------
% Neutral build for lookup
% -----------------------------
opts0 = struct();
opts0.supercell_size = [2 2 1];
opts0.verbose = true;
opts0.removeMolIDs = [];
opts0.activeMolIDs = [];

params0 = util.default_params();
result0 = calc.run_polarization_calc(crystal, model, opts0, params0);
sys0 = result0.sys;

molA = builder.find_unique_molecule_id(sys0, 1, [0 0 0]);
molB = builder.find_unique_molecule_id(sys0, 1, [1 0 0]);

opts = struct();
opts.supercell_size = [2 2 1];
opts.verbose = true;
opts.removeMolIDs = [];
opts.activeMolIDs = [molA, molB];

% -----------------------------
% Direct periodic solve
% -----------------------------
params_dir = util.default_params();
params_dir.ewald.mode = 'periodic_triclinic';
params_dir.ewald.auto = true;
params_dir.ewald.tol = 1e-8;
params_dir.ewald.boundary = 'tinfoil';
params_dir.ewald.rcut_fraction = 0.9;
params_dir.scf.solver = 'direct';
params_dir.scf.use_thole = false;

res_dir = calc.run_polarization_calc(crystal, model, opts, params_dir);

% -----------------------------
% Matrix-iterative periodic solve
% -----------------------------
params_iter = util.default_params();
params_iter.ewald.mode = 'periodic_triclinic';
params_iter.ewald.auto = true;
params_iter.ewald.tol = 1e-8;
params_iter.ewald.boundary = 'tinfoil';
params_iter.ewald.rcut_fraction = 0.9;
params_iter.scf.solver = 'matrix_iterative';
params_iter.scf.tol = 1e-12;
params_iter.scf.maxIter = 1000;
params_iter.scf.mixing = 0.35;
params_iter.scf.use_thole = false;

res_iter = calc.run_polarization_calc(crystal, model, opts, params_iter);

% -----------------------------
% Compare
% -----------------------------
dmu = res_iter.mu - res_dir.mu;
err = max(sqrt(sum(dmu.^2, 2)));

fprintf('Max |delta mu|         = %.16e\n', err);
fprintf('Direct Epol            = %.16e\n', res_dir.energy.polarization);
fprintf('Iterative Epol         = %.16e\n', res_iter.energy.polarization);
fprintf('Iterative converged    = %d\n', res_iter.scf.converged);
fprintf('Iterative iterations   = %d\n', res_iter.scf.nIter);

if isfield(res_dir, 'direct') && isfield(res_dir.direct, 'residual_norm')
    fprintf('Direct residual norm   = %.16e\n', res_dir.direct.residual_norm);
end