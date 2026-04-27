function test_io_site_classes()
%TEST_IO_SITE_CLASSES Verify graph-based site-class inference.
%
% We use a propene-like connectivity pattern:
%
%   C1H2 = C2H - C3H3
%
% The bond graph does not encode bond order, but the degrees are enough:
%
%   C1 degree 3 -> C_deg3
%   C2 degree 3 -> C_deg3
%   C3 degree 4 -> C_deg4
%
% Hydrogens inherit H_on_C_deg3 or H_on_C_deg4 from their bonded carbon.

species = {
    'C'  % 1 C1, alkene CH2 carbon, degree 3
    'C'  % 2 C2, alkene CH carbon, degree 3
    'C'  % 3 C3, methyl carbon, degree 4
    'H'  % 4 H on C1
    'H'  % 5 H on C1
    'H'  % 6 H on C2
    'H'  % 7 H on C3
    'H'  % 8 H on C3
    'H'  % 9 H on C3
};

n = numel(species);
A = false(n, n);

% Carbon skeleton: C1-C2-C3
A = local_add_bond(A, 1, 2);
A = local_add_bond(A, 2, 3);

% Hydrogens
A = local_add_bond(A, 1, 4);
A = local_add_bond(A, 1, 5);
A = local_add_bond(A, 2, 6);
A = local_add_bond(A, 3, 7);
A = local_add_bond(A, 3, 8);
A = local_add_bond(A, 3, 9);

site_class = io.infer_simple_site_classes(species, sparse(A));

expected = {
    'C_deg3'
    'C_deg3'
    'C_deg4'
    'H_on_C_deg3'
    'H_on_C_deg3'
    'H_on_C_deg3'
    'H_on_C_deg4'
    'H_on_C_deg4'
    'H_on_C_deg4'
};

assert(isequal(site_class, expected), ...
    'Propene-like site classes were inferred incorrectly.');

% Also check string input.
site_class_string = io.infer_simple_site_classes(string(species), sparse(A));

assert(isequal(site_class_string, expected), ...
    'String species input should produce the same site classes.');

% Unsupported carbon degree should fail.
A_bad = A;
A_bad(1, 4) = false;
A_bad(4, 1) = false;

didFail = false;
try
    io.infer_simple_site_classes(species, sparse(A_bad));
catch
    didFail = true;
end

assert(didFail, ...
    'Carbon degree 2 should fail for the simple site-class model.');

end

function A = local_add_bond(A, i, j)
A(i, j) = true;
A(j, i) = true;
end