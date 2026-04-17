%% run_test_import_supercell
clear; clc;

rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
filename   = fullfile(rootFolder, 'a_0.0_CONTCAR.vasp');

%% ------------------------------------------------------------------------
% Import crystal through the current io stack
crystal = io.import_contcar_as_crystal(filename, ...
    'BondScale', 1.20, ...
    'WrapFractional', true, ...
    'SortMolecules', false);

fprintf('Imported crystal:\n');
fprintf('  nSites   = %d\n', size(crystal.cart_coords, 1));
fprintf('  nMols    = %d\n', numel(unique(crystal.mol_id)));
fprintf('  lattice:\n');
disp(crystal.lattice)

%% ------------------------------------------------------------------------
% Build minimal unit struct expected by geom.build_supercell
unit = struct();
unit.lattice = crystal.lattice;
unit.frac_coords = crystal.frac_coords;
unit.cart_coords = crystal.cart_coords;
unit.mol_id = crystal.mol_id;
unit.site_label = crystal.site_label;
unit.site_type = crystal.site_type;

%% ------------------------------------------------------------------------
% Replicate
supercell_size = [2 2 1];
super = geom.build_supercell(unit, supercell_size);

fprintf('\nBuilt supercell [%d %d %d]:\n', supercell_size(1), supercell_size(2), supercell_size(3));
fprintf('  n_unit_sites = %d\n', super.n_unit_sites);
fprintf('  n_cells      = %d\n', super.n_cells);
fprintf('  n_sites      = %d\n', super.n_sites);

expected_n_sites = super.n_unit_sites * prod(supercell_size);
fprintf('  expected     = %d\n', expected_n_sites);

if super.n_sites ~= expected_n_sites
    error('Supercell site count mismatch.');
end

%% ------------------------------------------------------------------------
% Basic sanity checks
if size(super.frac_coords,1) ~= super.n_sites
    error('super.frac_coords row count mismatch.');
end

if size(super.cart_coords,1) ~= super.n_sites
    error('super.cart_coords row count mismatch.');
end

if size(super.cell_shift,1) ~= super.n_sites
    error('super.cell_shift row count mismatch.');
end

if numel(super.image_id) ~= super.n_sites
    error('super.image_id size mismatch.');
end

if numel(super.unit_site_index) ~= super.n_sites
    error('super.unit_site_index size mismatch.');
end

fprintf('\nUnique cell shifts:\n');
disp(unique(super.cell_shift, 'rows'))

fprintf('Unique inherited mol_ids in supercell = %d\n', numel(unique(super.mol_id)));

%% ------------------------------------------------------------------------
% Check coordinate consistency on a few sites
cart_from_frac = geom.frac_to_cart(super.frac_coords, super.lattice);
err = max(abs(cart_from_frac(:) - super.cart_coords(:)));
fprintf('Max |frac_to_cart(super.frac_coords) - super.cart_coords| = %.3e\n', err);

if err > 1e-10
    warning('Supercell Cartesian/fractional coordinates are not perfectly consistent.');
end

%% ------------------------------------------------------------------------
% Show first few rows for inspection
nShow = min(10, super.n_sites);
T = table( ...
    (1:nShow).', ...
    super.unit_site_index(1:nShow), ...
    super.image_id(1:nShow), ...
    super.cell_shift(1:nShow,1), ...
    super.cell_shift(1:nShow,2), ...
    super.cell_shift(1:nShow,3), ...
    super.mol_id(1:nShow), ...
    string(super.site_type(1:nShow)), ...
    string(super.site_label(1:nShow)), ...
    'VariableNames', { ...
        'site', 'unit_site_index', 'image_id', ...
        'ix', 'iy', 'iz', 'mol_id', 'site_type', 'site_label'});

disp(T)

fprintf('\nCheckpoint passed.\n');