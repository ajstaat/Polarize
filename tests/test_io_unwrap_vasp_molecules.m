function test_io_unwrap_vasp_molecules()
%TEST_IO_UNWRAP_VASP_MOLECULES Verify PBC molecule unwrapping.

H = [
    10.0 0.0 0.0
     0.0 10.0 0.0
     0.0 0.0 10.0
];

S = struct();
S.comment = 'Tiny boundary-crossing molecule';
S.scale = 1.0;
S.lattice = H;
S.species = {'C'; 'C'; 'H'};
S.frac = [
    0.95 0.50 0.50
    0.05 0.50 0.50
    0.50 0.50 0.50
];
S.cart = S.frac * H;
S.coord_type = 'direct';
S.natoms = numel(S.species);

[A, ~] = io.build_pbc_bond_graph(S.species, S.frac, S.lattice, 1.20, ...
    'Method', 'cell_list');

% The two carbons should form one connected component crossing the boundary.
mol = io.unwrap_vasp_molecule(S, [1; 2], A);

assert(isequal(mol.indices, [1; 2]), ...
    'Unwrapped molecule should preserve requested atom indices.');

assert(numel(mol.labels) == 2, ...
    'Unwrapped molecule should contain two labels.');

assert(isequal(size(mol.fracWrapped), [2 3]), ...
    'fracWrapped should be 2 x 3.');

assert(isequal(size(mol.fracUnwrapped), [2 3]), ...
    'fracUnwrapped should be 2 x 3.');

assert(isequal(size(mol.cart), [2 3]), ...
    'cart should be 2 x 3.');

d12 = norm(mol.cart(2,:) - mol.cart(1,:));

assert(abs(d12 - 1.0) < 1e-12, ...
    'Unwrapped boundary-crossing molecule should have 1.0 Cartesian separation.');

% The unwrapped fractional displacement from atom 1 to atom 2 should be
% +0.1 along x, not -0.9.
df12 = mol.fracUnwrapped(2,:) - mol.fracUnwrapped(1,:);

assert(norm(df12 - [0.10 0.00 0.00], 'fro') < 1e-12, ...
    'Unwrapped fractional displacement should use the minimum image.');

% Now test unwrap_all with a supplied bond graph.
molecules = io.unwrap_all_contcar_molecules(S, ...
    'BondScale', 1.20, ...
    'BondGraph', A, ...
    'SortMolecules', false);

assert(numel(molecules) == 2, ...
    'Expected two connected components: C-C molecule and isolated H.');

sizes = cellfun(@(m) numel(m.indices), molecules);

assert(any(sizes == 2), ...
    'One molecule/component should contain the bonded C-C pair.');

assert(any(sizes == 1), ...
    'One molecule/component should contain the isolated H.');

idx2 = find(sizes == 2, 1, 'first');
mol2 = molecules{idx2};

d = norm(mol2.cart(2,:) - mol2.cart(1,:));

assert(abs(d - 1.0) < 1e-12, ...
    'unwrap_all_contcar_molecules should unwrap the C-C pair coherently.');

% Also test that unwrap_all can build the bond graph internally.
molecules2 = io.unwrap_all_contcar_molecules(S, ...
    'BondScale', 1.20, ...
    'SortMolecules', false);

sizes2 = sort(cellfun(@(m) numel(m.indices), molecules2));

assert(isequal(sizes2(:), [1; 2]), ...
    'Internal bond-graph path should produce one 1-site and one 2-site component.');

end