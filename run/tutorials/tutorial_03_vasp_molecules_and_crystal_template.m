%% tutorial_03_vasp_molecules_and_crystal_template
% This tutorial demonstrates the molecule-aware VASP IO layer in Polarize:
%
%   1. reading a small VASP POSCAR-style file
%   2. building a PBC-aware covalent bond graph
%   3. unwrapping connected molecules under periodic boundary conditions
%   4. importing the same file as a crystal template for the builder layer
%
% The example molecule is propene-like:
%
%   C1H2 = C2H - C3H3
%
% The bond graph does not encode bond order, but the local degrees are:
%
%   C1 degree 3 -> C_deg3
%   C2 degree 3 -> C_deg3
%   C3 degree 4 -> C_deg4
%
% Hydrogens are classified according to the carbon they are bonded to:
%
%   H_on_C_deg3
%   H_on_C_deg4
%
% This tutorial assumes the MATLAB project has already added src/ to the path.

clear; clc;

fprintf('\n============================================================\n');
fprintf('Tutorial 03: VASP molecules and crystal template\n');
fprintf('============================================================\n');

%% 1. Write a tiny temporary POSCAR-style file

tmpDir = tempname;
mkdir(tmpDir);
cleanup = onCleanup(@() local_cleanup(tmpDir));

filename = fullfile(tmpDir, 'POSCAR_propene_like');

% Coordinates are Cartesian Angstrom in a large orthorhombic box.
%
% Atom order:
%   1 C1
%   2 C2
%   3 C3
%   4 H on C1
%   5 H on C1
%   6 H on C2
%   7 H on C3
%   8 H on C3
%   9 H on C3

xyz = [
    5.00  5.00  5.00   % C1
    6.45  5.00  5.00   % C2
    7.90  5.00  5.00   % C3
    4.55  5.95  5.00   % H on C1
    4.55  4.05  5.00   % H on C1
    6.45  6.09  5.00   % H on C2
    8.35  5.95  5.00   % H on C3
    8.35  4.05  5.00   % H on C3
    7.90  5.00  6.09   % H on C3
];

fid = fopen(filename, 'w');
if fid < 0
    error('Failed to open temporary POSCAR file.');
end

fprintf(fid, 'Propene-like tutorial POSCAR\n');
fprintf(fid, '1.0\n');
fprintf(fid, '20.0 0.0 0.0\n');
fprintf(fid, '0.0 20.0 0.0\n');
fprintf(fid, '0.0 0.0 20.0\n');
fprintf(fid, 'C H\n');
fprintf(fid, '3 6\n');
fprintf(fid, 'Cartesian\n');

for i = 1:size(xyz, 1)
    fprintf(fid, '%.10f %.10f %.10f\n', xyz(i,1), xyz(i,2), xyz(i,3));
end

fclose(fid);

fprintf('\nWrote temporary POSCAR:\n  %s\n', filename);

%% 2. Read the raw VASP structure

S = io.read_vasp_structure(filename);

fprintf('\nRaw VASP parser output:\n');
fprintf('  comment    = %s\n', S.comment);
fprintf('  coord_type = %s\n', S.coord_type);
fprintf('  natoms     = %d\n', S.natoms);
fprintf('  lattice    = [%d x %d]\n', size(S.lattice,1), size(S.lattice,2));
fprintf('  frac       = [%d x %d]\n', size(S.frac,1), size(S.frac,2));
fprintf('  cart       = [%d x %d]\n', size(S.cart,1), size(S.cart,2));

fprintf('\nSpecies labels:\n');
disp(S.species(:).');

%% 3. Build the PBC-aware covalent bond graph

bondScale = 1.20;

[A, Dcart] = io.build_pbc_bond_graph( ...
    S.species, S.frac, S.lattice, bondScale, ...
    'ReturnDcart', true, ...
    'Method', 'cell_list');

nBonds = nnz(triu(A, 1));

fprintf('\nPBC bond graph:\n');
fprintf('  bond scale     = %.2f\n', bondScale);
fprintf('  accepted bonds = %d\n', nBonds);

fprintf('\nBond list:\n');
fprintf('  i  j   type      distance / Angstrom\n');
fprintf('  -- --  --------  -------------------\n');

[ii, jj] = find(triu(A, 1));

for k = 1:numel(ii)
    i = ii(k);
    j = jj(k);

    pairLabel = sprintf('%s-%s', S.species{i}, S.species{j});

    fprintf('  %2d %2d  %-8s  %10.4f\n', ...
        i, j, pairLabel, Dcart(i,j));
end

%% 4. Unwrap all connected molecular components

molecules = io.unwrap_all_contcar_molecules(S, ...
    'BondScale', bondScale, ...
    'BondGraph', A, ...
    'SortMolecules', false);

fprintf('\nUnwrapped molecule components:\n');
fprintf('  number of components = %d\n', numel(molecules));

for m = 1:numel(molecules)
    mol = molecules{m};

    fprintf('\n  Molecule/component %d:\n', m);
    fprintf('    n sites = %d\n', numel(mol.indices));
    fprintf('    indices = ');
    fprintf('%d ', mol.indices);
    fprintf('\n');

    fprintf('    labels  = ');
    fprintf('%s ', mol.labels{:});
    fprintf('\n');

    fprintf('    first unwrapped Cartesian coordinate / Angstrom:\n');
    disp(mol.cart(1,:));
end

%% 5. Import as a builder-ready crystal template

crystal = io.import_contcar_as_crystal(filename, ...
    'BondScale', bondScale, ...
    'SortMolecules', false);

fprintf('\nCrystal template:\n');
fprintf('  nSites     = %d\n', crystal.nSites);
fprintf('  nBaseMols  = %d\n', crystal.nBaseMols);
fprintf('  length unit = %s\n', crystal.units.length);

fprintf('\nCrystal fields used by builder layer:\n');
fprintf('  lattice      = [%d x %d]  bohr\n', ...
    size(crystal.lattice,1), size(crystal.lattice,2));
fprintf('  frac_coords  = [%d x %d]  dimensionless\n', ...
    size(crystal.frac_coords,1), size(crystal.frac_coords,2));
fprintf('  cart_coords  = [%d x %d]  bohr\n', ...
    size(crystal.cart_coords,1), size(crystal.cart_coords,2));
fprintf('  base_mol_id  = [%d x %d]\n', ...
    size(crystal.base_mol_id,1), size(crystal.base_mol_id,2));
fprintf('  site_type    = [%d x %d]\n', ...
    size(crystal.site_type,1), size(crystal.site_type,2));
fprintf('  site_class   = [%d x %d]\n', ...
    size(crystal.site_class,1), size(crystal.site_class,2));
fprintf('  site_label   = [%d x %d]\n', ...
    size(crystal.site_label,1), size(crystal.site_label,2));

fprintf('\nPer-site crystal metadata:\n');
fprintf('  idx  type  base_mol  label  class\n');
fprintf('  ---  ----  --------  -----  --------------\n');

for i = 1:crystal.nSites
    fprintf('  %3d  %-4s  %8d  %-5s  %s\n', ...
        i, ...
        crystal.site_type{i}, ...
        crystal.base_mol_id(i), ...
        crystal.site_label{i}, ...
        crystal.site_class{i});
end

fprintf('\nCoordinate unit conversion example:\n');
fprintf('  raw VASP Cartesian coordinate for site 2 / Angstrom:\n');
disp(S.cart(2,:));

fprintf('  crystal Cartesian coordinate for site 2 / bohr:\n');
disp(crystal.cart_coords(2,:));

fprintf('\nTutorial 03 completed successfully.\n');

function local_cleanup(tmpDir)
if exist(tmpDir, 'dir')
    rmdir(tmpDir, 's');
end
end