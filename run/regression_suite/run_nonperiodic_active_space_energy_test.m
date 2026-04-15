function run_nonperiodic_active_space_energy_test
clc;

fprintf('=== nonperiodic active-space energy/decomposition test ===\n');

% -------------------------------------------------------------------------
% Small case patterned after your run scripts
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

params = util.default_params();
params.scf.use_thole = true;
params.scf.tol = 1e-12;
params.scf.maxIter = 1000;
params.scf.mixing = 0.35;
params.scf.omega = 1.0;
params.output.computeEnergy = true;
params.output.computeEnergyByMolecule = true;
params.output.computeDipoleDipoleDecomp = true;

params_iter = params;
params_iter.scf.solver = 'matrix_iterative';

params_dir = params;
params_dir.scf.solver = 'direct';

params_sor = params;
params_sor.scf.solver = 'sor';

res_iter = calc.run_polarization_calc(crystal, model, opts, params_iter);
res_dir  = calc.run_polarization_calc(crystal, model, opts, params_dir);
res_sor  = calc.run_polarization_calc(crystal, model, opts, params_sor);

% -------------------------------------------------------------------------
% 1) Dipole agreement
% -------------------------------------------------------------------------
fprintf('\n--- dipole agreement ---\n');

err_iter_dir = max(sqrt(sum((res_iter.mu - res_dir.mu).^2, 2)));
err_sor_dir  = max(sqrt(sum((res_sor.mu  - res_dir.mu).^2, 2)));

fprintf('max |mu_iter - mu_direct| = %.16e\n', err_iter_dir);
fprintf('max |mu_sor  - mu_direct| = %.16e\n', err_sor_dir);

assert(err_iter_dir < 1e-8, 'matrix_iterative mu disagrees with direct');
assert(err_sor_dir  < 1e-8, 'sor mu disagrees with direct');

% -------------------------------------------------------------------------
% 2) Total energy agreement across solvers
% -------------------------------------------------------------------------
fprintf('\n--- total energy agreement ---\n');

assert(isfield(res_iter, 'energy') && isfield(res_iter.energy, 'total'), 'iter energy.total missing');
assert(isfield(res_dir,  'energy') && isfield(res_dir.energy,  'total'), 'dir energy.total missing');
assert(isfield(res_sor,  'energy') && isfield(res_sor.energy,  'total'), 'sor energy.total missing');

E_iter = res_iter.energy.total;
E_dir  = res_dir.energy.total;
E_sor  = res_sor.energy.total;

fprintf('E_total(iter) = %.16e\n', E_iter);
fprintf('E_total(dir ) = %.16e\n', E_dir);
fprintf('E_total(sor ) = %.16e\n', E_sor);

fprintf('|E_iter - E_dir| = %.16e\n', abs(E_iter - E_dir));
fprintf('|E_sor  - E_dir| = %.16e\n', abs(E_sor - E_dir));

assert(abs(E_iter - E_dir) < 1e-10, 'matrix_iterative total energy disagrees with direct');
assert(abs(E_sor  - E_dir) < 1e-10, 'sor total energy disagrees with direct');

% -------------------------------------------------------------------------
% 3) Per-molecule energy table agreement
% -------------------------------------------------------------------------
fprintf('\n--- per-molecule energy agreement ---\n');

tbl_iter = res_iter.energy_by_molecule;
tbl_dir  = res_dir.energy_by_molecule;
tbl_sor  = res_sor.energy_by_molecule;

assert(istable(tbl_iter), 'iter energy_by_molecule is not a table');
assert(istable(tbl_dir),  'dir energy_by_molecule is not a table');
assert(istable(tbl_sor),  'sor energy_by_molecule is not a table');

assert(height(tbl_iter) == height(tbl_dir), 'iter/dir molecule table height mismatch');
assert(height(tbl_sor)  == height(tbl_dir), 'sor/dir molecule table height mismatch');

if any(strcmp(tbl_iter.Properties.VariableNames, 'mol_id'))
    tbl_iter = sortrows(tbl_iter, 'mol_id');
    tbl_dir  = sortrows(tbl_dir,  'mol_id');
    tbl_sor  = sortrows(tbl_sor,  'mol_id');
end

vars_iter = tbl_iter.Properties.VariableNames;
vars_dir  = tbl_dir.Properties.VariableNames;
vars_sor  = tbl_sor.Properties.VariableNames;

assert(isequal(vars_iter, vars_dir), 'iter/dir molecule table variable names mismatch');
assert(isequal(vars_sor,  vars_dir), 'sor/dir molecule table variable names mismatch');

for c = 1:numel(vars_dir)
    name = vars_dir{c};

    if isnumeric(tbl_dir.(name))
        diff_iter = max(abs(tbl_iter.(name) - tbl_dir.(name)));
        diff_sor  = max(abs(tbl_sor.(name)  - tbl_dir.(name)));

        fprintf('column %-24s | max|iter-dir| = %.3e | max|sor-dir| = %.3e\n', ...
            name, diff_iter, diff_sor);

        assert(diff_iter < 1e-10, 'per-molecule energy mismatch: iter vs dir');
        assert(diff_sor  < 1e-10, 'per-molecule energy mismatch: sor vs dir');
    else
        assert(isequal(tbl_iter.(name), tbl_dir.(name)), 'nonnumeric per-molecule column mismatch: iter vs dir');
        assert(isequal(tbl_sor.(name),  tbl_dir.(name)), 'nonnumeric per-molecule column mismatch: sor vs dir');
    end
end

if any(strcmp(tbl_dir.Properties.VariableNames, 'total'))
    sum_mol_dir = sum(tbl_dir.total);
    fprintf('sum(molecule total)    = %.16e\n', sum_mol_dir);
    fprintf('global total energy    = %.16e\n', E_dir);
    fprintf('|sum_mol - global|     = %.16e\n', abs(sum_mol_dir - E_dir));
end

% -------------------------------------------------------------------------
% 4) Dipole-dipole decomposition agreement
% -------------------------------------------------------------------------
fprintf('\n--- dipole-dipole decomposition agreement ---\n');

dd_iter = res_iter.dipole_dipole_decomposition;
dd_dir  = res_dir.dipole_dipole_decomposition;
dd_sor  = res_sor.dipole_dipole_decomposition;

assert(isstruct(dd_iter), 'iter dipole_dipole_decomposition is not a struct');
assert(isstruct(dd_dir),  'dir dipole_dipole_decomposition is not a struct');
assert(isstruct(dd_sor),  'sor dipole_dipole_decomposition is not a struct');

if isfield(dd_iter, 'total') && isfield(dd_dir, 'total') && isfield(dd_sor, 'total')
    fprintf('dd total(iter) = %.16e\n', dd_iter.total);
    fprintf('dd total(dir ) = %.16e\n', dd_dir.total);
    fprintf('dd total(sor ) = %.16e\n', dd_sor.total);

    fprintf('|dd_iter - dd_dir| = %.16e\n', abs(dd_iter.total - dd_dir.total));
    fprintf('|dd_sor  - dd_dir| = %.16e\n', abs(dd_sor.total  - dd_dir.total));

    assert(abs(dd_iter.total - dd_dir.total) < 1e-10, 'dd total mismatch: iter vs dir');
    assert(abs(dd_sor.total  - dd_dir.total) < 1e-10, 'dd total mismatch: sor vs dir');
end

has_pair_arrays = isfield(dd_dir, 'energy') && isfield(dd_iter, 'energy') && isfield(dd_sor, 'energy');
if has_pair_arrays
    assert(numel(dd_iter.energy) == numel(dd_dir.energy), 'iter/dir dd pair count mismatch');
    assert(numel(dd_sor.energy)  == numel(dd_dir.energy), 'sor/dir dd pair count mismatch');

    if isfield(dd_dir, 'site_i') && isfield(dd_dir, 'site_j')
        key_dir  = [dd_dir.site_i(:),  dd_dir.site_j(:)];
        key_iter = [dd_iter.site_i(:), dd_iter.site_j(:)];
        key_sor  = [dd_sor.site_i(:),  dd_sor.site_j(:)];

        [key_dir, ord_dir] = sortrows(key_dir);
        [key_iter, ord_iter] = sortrows(key_iter);
        [key_sor, ord_sor] = sortrows(key_sor);

        assert(isequal(key_iter, key_dir), 'iter dd pair keys mismatch');
        assert(isequal(key_sor,  key_dir), 'sor dd pair keys mismatch');

        e_dir  = dd_dir.energy(ord_dir);
        e_iter = dd_iter.energy(ord_iter);
        e_sor  = dd_sor.energy(ord_sor);
    else
        e_dir  = dd_dir.energy(:);
        e_iter = dd_iter.energy(:);
        e_sor  = dd_sor.energy(:);
    end

    fprintf('max |dd_iter - dd_dir| = %.16e\n', max(abs(e_iter - e_dir)));
    fprintf('max |dd_sor  - dd_dir| = %.16e\n', max(abs(e_sor  - e_dir)));

    assert(max(abs(e_iter - e_dir)) < 1e-10, 'pairwise dd mismatch: iter vs dir');
    assert(max(abs(e_sor  - e_dir)) < 1e-10, 'pairwise dd mismatch: sor vs dir');

    if isfield(dd_dir, 'total')
        fprintf('sum(dd pair energies)  = %.16e\n', sum(e_dir));
        fprintf('dd total              = %.16e\n', dd_dir.total);
        fprintf('|sum_pairs - total|   = %.16e\n', abs(sum(e_dir) - dd_dir.total));
    end
end

% -------------------------------------------------------------------------
% 5) Optional direct active-space consistency checks
% -------------------------------------------------------------------------
fprintf('\n--- active-space internal consistency ---\n');

sys = builder.make_crystal_system(crystal, model, opts);
sys = builder.assign_point_charges(sys, model.charge_pattern);
Eext = calc.compute_external_field(sys, params);
problem = thole.prepare_scf_problem(sys, Eext, params.scf);
[Tpol, ~] = ewald.assemble_nonperiodic_interaction_matrix(sys, problem, params.ewald, params.scf);

mu_dir = res_dir.mu;

if exist('calc.compute_total_energy_active_space', 'file') == 2
    e_active = calc.compute_total_energy_active_space(sys, problem, mu_dir, Eext, Tpol);
    fprintf('direct active-space total = %.16e\n', e_active.total);
    fprintf('driver total              = %.16e\n', E_dir);
    fprintf('|active - driver|         = %.16e\n', abs(e_active.total - E_dir));
    assert(abs(e_active.total - E_dir) < 1e-10, 'active-space total-energy helper disagrees with driver');
end

if exist('calc.compute_total_energy_by_molecule_active_space', 'file') == 2
    tbl_active = calc.compute_total_energy_by_molecule_active_space(sys, problem, mu_dir, Eext, Tpol);
    if any(strcmp(tbl_active.Properties.VariableNames, 'mol_id'))
        tbl_active = sortrows(tbl_active, 'mol_id');
        tbl_dir2 = sortrows(tbl_dir, 'mol_id');
    else
        tbl_dir2 = tbl_dir;
    end

    for c = 1:numel(tbl_dir2.Properties.VariableNames)
        name = tbl_dir2.Properties.VariableNames{c};
        if isnumeric(tbl_dir2.(name))
            assert(max(abs(tbl_active.(name) - tbl_dir2.(name))) < 1e-10, ...
                'active-space per-molecule helper disagrees with driver');
        end
    end
end

fprintf('\nPASS: active-space energy/decomposition test completed.\n');
end