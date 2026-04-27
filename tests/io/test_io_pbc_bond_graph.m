function test_io_pbc_bond_graph()
%TEST_IO_PBC_BOND_GRAPH Verify PBC-aware covalent bond graph construction.
%
% This test checks that io.build_pbc_bond_graph uses the shared spatial
% index with the project row-lattice convention:
%
%   cart = frac * H

H = [
    10.0 0.0 0.0
     0.0 10.0 0.0
     0.0 0.0 10.0
];

% Two carbons across the x boundary:
% wrapped separation is 9.0, minimum-image separation is 1.0.
%
% Third atom is distant from both carbons.
species = {'C'; 'C'; 'H'};

frac = [
    0.95 0.50 0.50
    0.05 0.50 0.50
    0.50 0.50 0.50
];

[A, D] = io.build_pbc_bond_graph(species, frac, H, 1.20, ...
    'ReturnDcart', true, ...
    'Method', 'cell_list');

assert(isequal(size(A), [3 3]), ...
    'Adjacency matrix should be N x N.');

assert(issparse(A), ...
    'Adjacency matrix should be sparse.');

assert(isequal(A, A.'), ...
    'Adjacency matrix should be symmetric.');

assert(~any(diag(A)), ...
    'Adjacency matrix diagonal should be false.');

assert(A(1,2), ...
    'Boundary-crossing C-C pair should be bonded by minimum image.');

assert(abs(D(1,2) - 1.0) < 1e-12, ...
    'Minimum-image distance across boundary should be 1.0.');

assert(~A(1,3) && ~A(2,3), ...
    'Distant C-H pairs should not be bonded in this test.');

% Compare against the brute-force spatial-index backend. This does not
% replace the performant backend; it verifies backend consistency.
[A_bf, D_bf] = io.build_pbc_bond_graph(species, frac, H, 1.20, ...
    'ReturnDcart', true, ...
    'Method', 'bruteforce');

assert(isequal(A, A_bf), ...
    'cell_list and bruteforce bond graphs should match.');

assert(norm(D - D_bf, 'fro') < 1e-12, ...
    'cell_list and bruteforce distance matrices should match.');

% H-H exclusion check.
speciesHH = {'H'; 'H'};

fracHH = [
    0.00 0.00 0.00
    0.05 0.00 0.00
];

[A_HH, ~] = io.build_pbc_bond_graph(speciesHH, fracHH, H, 1.20, ...
    'Method', 'cell_list');

assert(~A_HH(1,2), ...
    'H-H pairs should be ignored by bond_graph_tools.ignore_pair.');

% Sheared row-lattice case. Cubic cells can hide row/column mistakes, so
% this explicitly tests a non-orthogonal row-vector lattice.
Hshear = [
    10.0  0.0  0.0
     2.0 10.0  0.0
     0.0  0.0 10.0
];

[A_shear, D_shear] = io.build_pbc_bond_graph(species, frac, Hshear, 1.20, ...
    'ReturnDcart', true, ...
    'Method', 'cell_list');

assert(A_shear(1,2), ...
    'Boundary-crossing C-C pair should be bonded in sheared cell.');

assert(abs(D_shear(1,2) - 1.0) < 1e-12, ...
    'Sheared-cell minimum-image distance should be 1.0.');

assert(~A_shear(1,3) && ~A_shear(2,3), ...
    'Distant C-H pairs should not be bonded in sheared-cell test.');

end