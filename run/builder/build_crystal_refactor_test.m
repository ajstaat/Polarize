clear; clc; close all;
fprintf('=== build crystal system smoke test ===\n');

cfg = struct();
cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');
cfg.supercellSize = [5 11 3];
cfg.bondScale = 1.20;
cfg.verbose = true;

fprintf('\n[setup] importing crystal template...\n');
crystal = io.import_contcar_as_crystal(cfg.filename, ...
    'BondScale', cfg.bondScale, ...
    'WrapFractional', true, ...
    'SortMolecules', false);

fprintf('Imported crystal template:\n');
fprintf('  nSites    = %d\n', size(crystal.cart_coords, 1));
fprintf('  nBaseMols = %d\n', numel(unique(crystal.base_mol_id)));

model = struct();
model.polarizable_classes = { ...
    'H_on_C_deg3', ...
    'H_on_C_deg4', ...
    'C_deg3', ...
    'C_deg4', ...
    'N', ...
    'O'};
model.alpha_by_class = struct( ...
    'H_on_C_deg3', 0.496, ...
    'H_on_C_deg4', 0.696, ...
    'C_deg3', 1.334, ...
    'C_deg4', 1.750, ...
    'N', 1.073, ...
    'O', 0.837);
model.thole_a = 0.39;

fprintf('\n[setup] building crystal system...\n');
buildOpts = struct();
buildOpts.supercell_size = cfg.supercellSize;
buildOpts.bondScale = cfg.bondScale;
buildOpts.verbose = cfg.verbose;
buildOpts.removeMolIDs = [];
buildOpts.activeMolIDs = [];

sys = builder.make_crystal_system(crystal, model, buildOpts);

fprintf('\nBuilt working system:\n');
fprintf('  n_unit_sites = %d\n', sys.n_unit_sites);
fprintf('  n_cells      = %d\n', sys.n_cells);
fprintf('  n_sites      = %d\n', sys.n_sites);

assert(sys.n_unit_sites > 0, 'n_unit_sites should be positive.');
assert(sys.n_cells == prod(cfg.supercellSize), 'n_cells mismatch.');
assert(sys.n_sites > 0, 'n_sites should be positive.');

if isfield(sys, 'molecule_id')
    molIDs = unique(sys.molecule_id(:));
    molIDs = molIDs(molIDs > 0);
    fprintf('  n_molecules   = %d\n', numel(molIDs));
    assert(~isempty(molIDs), 'No molecules identified in built system.');
end

fprintf('\nSmoke test passed.\n');