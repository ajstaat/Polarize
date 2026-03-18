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
res0 = calc.run_polarization_calc(crystal, model, opts0, params0);
sys0 = res0.sys;

molA = builder.find_unique_molecule_id(sys0, 1, [0 0 0]);
molB = builder.find_unique_molecule_id(sys0, 1, [1 0 0]);

opts = struct();
opts.supercell_size = [2 2 1];
opts.verbose = true;
opts.removeMolIDs = [];
opts.activeMolIDs = [molA, molB];

% -----------------------------
% Bare periodic operator
% -----------------------------
params_bare = util.default_params();
params_bare.ewald.mode = 'periodic_triclinic';
params_bare.ewald.auto = true;
params_bare.ewald.tol = 1e-8;
params_bare.ewald.boundary = 'tinfoil';
params_bare.ewald.rcut_fraction = 0.9;
params_bare.ewald.use_thole_real_space = false;
params_bare.scf.solver = 'direct';

res_bare = calc.run_polarization_calc(crystal, model, opts, params_bare);

% -----------------------------
% Periodic operator with real-space Thole correction
% -----------------------------
params_thole = util.default_params();
params_thole.ewald.mode = 'periodic_triclinic';
params_thole.ewald.auto = true;
params_thole.ewald.tol = 1e-8;
params_thole.ewald.boundary = 'tinfoil';
params_thole.ewald.rcut_fraction = 0.9;
params_thole.ewald.use_thole_real_space = true;
params_thole.ewald.thole_a = model.thole_a;
params_thole.scf.solver = 'direct';

res_thole = calc.run_polarization_calc(crystal, model, opts, params_thole);

fprintf('Bare periodic total energy         = %+ .15e\n', res_bare.energy.total);
fprintf('Thole-corrected periodic total     = %+ .15e\n', res_thole.energy.total);
fprintf('|delta U|                          = %.3e\n', ...
    abs(res_thole.energy.total - res_bare.energy.total));

fprintf('\nBare periodic dipole-dipole decomposition:\n');
disp(res_bare.dipole_dipole_decomposition)

fprintf('\nThole-corrected periodic dipole-dipole decomposition:\n');
disp(res_thole.dipole_dipole_decomposition)

fprintf('\nBare periodic induced dipoles:\n');
disp(res_bare.mu)

fprintf('\nThole-corrected periodic induced dipoles:\n');
disp(res_thole.mu)