%% run_test_io_import
clear; clc;

rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
filename   = fullfile(rootFolder, 'a_0.0_CONTCAR.vasp');

S = io.read_vasp_structure(filename);
disp(S.comment)
fprintf('natoms = %d\n', S.natoms);

mols = io.unwrap_all_contcar_molecules(S, ...
    'BondScale', 1.20, ...
    'SortMolecules', false);

fprintf('nmolecules = %d\n', numel(mols));
fprintf('atoms in first molecule = %d\n', numel(mols{1}.indices));

crystal = io.import_contcar_as_crystal(filename, ...
    'BondScale', 1.20, ...
    'WrapFractional', true, ...
    'SortMolecules', false);

fprintf('crystal sites = %d\n', size(crystal.cart_coords,1));
fprintf('unique mol IDs = %d\n', numel(unique(crystal.mol_id)));
disp(crystal.site_label(1:min(10,numel(crystal.site_label))))