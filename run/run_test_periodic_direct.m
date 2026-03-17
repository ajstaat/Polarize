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
% Build neutral system for molecule lookup
% -----------------------------
opts0 = struct();
opts0.supercell_size = [2 2 1];
opts0.verbose = true;
opts0.removeMolIDs = [];
opts0.activeMolIDs = [];

params0 = util.default_params();
result0 = calc.run_polarization_calc(crystal, model, opts0, params0);
sys0 = result0.sys;

Tmol = builder.list_molecules(sys0);
disp(Tmol)

molA = builder.find_unique_molecule_id(sys0, 1, [0 0 0]);
molB = builder.find_unique_molecule_id(sys0, 1, [1 0 0]);

% -----------------------------
% Final periodic run
% -----------------------------
opts = struct();
opts.supercell_size = [2 2 1];
opts.verbose = true;
opts.removeMolIDs = [];
opts.activeMolIDs = [molA, molB];

params = util.default_params();

% Important: use periodic operator
params.ewald.mode = 'periodic_triclinic';

% Safer first test: auto-select Ewald parameters from the actual supercell
params.ewald.auto = true;
params.ewald.tol = 1e-8;
params.ewald.boundary = 'tinfoil';
params.ewald.rcut_fraction = 0.9;

% Use direct solver first
params.scf.solver = 'direct';

% For now keep this off in the periodic operator path,
% since the current periodic T is not yet Thole-damped.
params.scf.use_thole = false;

result = calc.run_polarization_calc(crystal, model, opts, params);

disp(result.status)
disp(result.message)

disp('SCF converged:')
disp(result.scf.converged)

disp('Ewald metadata:')
disp(result.Tmeta)

disp('Site charges:')
disp(result.sys.site_charge(:).')

disp('External field:')
disp(result.Eext)

disp('Induced dipoles:')
disp(result.mu)

disp('Polarization energy:')
disp(result.energy.polarization)

if isfield(result, 'direct') && isfield(result.direct, 'residual_norm')
    fprintf('Direct residual norm = %.16e\n', result.direct.residual_norm);
end
if isfield(result, 'direct') && isfield(result.direct, 'condition_estimate')
    fprintf('Condition estimate   = %.16e\n', result.direct.condition_estimate);
end