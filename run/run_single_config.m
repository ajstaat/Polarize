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
% Neutral build for molecule lookup
% -----------------------------
opts0 = struct();
opts0.supercell_size = [2 2 1];
opts0.verbose = true;
opts0.removeMolIDs = [];
opts0.activeMolIDs = [];

params = util.default_params();
params.scf.solver = 'matrix_iterative';
params.scf.tol = 1e-10;
params.scf.maxIter = 500;
params.scf.mixing = 0.35;
params.scf.use_thole = true;

result0 = calc.run_polarization_calc(crystal, model, opts0, params);
sys0 = result0.sys;

Tmol = builder.list_molecules(sys0);
disp(Tmol)

molA = builder.find_unique_molecule_id(sys0, 1, [0 0 0]);
molB = builder.find_unique_molecule_id(sys0, 1, [1 0 0]);

% -----------------------------
% Final run
% -----------------------------
opts = struct();
opts.supercell_size = [2 2 1];
opts.verbose = true;
opts.removeMolIDs = [];
opts.activeMolIDs = [molA, molB];

result = calc.run_polarization_calc(crystal, model, opts, params);

disp(result.status)
disp(result.message)

disp('SCF converged:')
disp(result.scf.converged)

disp('Iterations:')
disp(result.scf.nIter)

disp('T size:')
disp(size(result.T))

disp('Induced dipoles:')
disp(result.mu)

disp('Polarization energy:')
disp(result.energy.polarization)