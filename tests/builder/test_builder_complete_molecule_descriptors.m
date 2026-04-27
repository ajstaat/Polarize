function test_builder_complete_molecule_descriptors()
%TEST_BUILDER_COMPLETE_MOLECULE_DESCRIPTORS Verify descriptor table values.

sys = local_make_descriptor_test_sys();

desc = builder.complete_molecule_descriptors_relative_to_reference( ...
    sys, 1, ...
    'StackAxis', [0 1 0], ...
    'IncludeReference', false, ...
    'IncludeNormals', true);

T = desc.table;

assert(desc.reference_mol_id == 1, ...
    'Unexpected reference molecule ID.');

assert(isequal(desc.stack_axis_hat, [0 1 0]), ...
    'Unexpected stack-axis unit vector.');

assert(height(T) == 4, ...
    'Expected four complete non-reference molecules in descriptor table.');

assert(~ismember(6, T.molecule_id), ...
    'Incomplete molecule should be excluded from descriptor table.');

% Molecule 2 is same-stack +b: dr = [0 2 0].
row2 = find(T.molecule_id == 2, 1, 'first');

assert(~isempty(row2), 'Molecule 2 missing from descriptor table.');
assert(abs(T.d_par(row2) - 2) < 1e-12, ...
    'Molecule 2 should have d_par = +2.');
assert(abs(T.d_perp(row2)) < 1e-12, ...
    'Molecule 2 should have d_perp = 0.');
assert(abs(T.distance(row2) - 2) < 1e-12, ...
    'Molecule 2 should have distance = 2.');

% Molecule 4 is side-stack +a: dr = [3 0 0].
row4 = find(T.molecule_id == 4, 1, 'first');

assert(~isempty(row4), 'Molecule 4 missing from descriptor table.');
assert(abs(T.d_par(row4)) < 1e-12, ...
    'Molecule 4 should have d_par = 0.');
assert(abs(T.d_perp(row4) - 3) < 1e-12, ...
    'Molecule 4 should have d_perp = 3.');
assert(norm([T.d_perp_x(row4), T.d_perp_y(row4), T.d_perp_z(row4)] - [3 0 0]) < 1e-12, ...
    'Molecule 4 should have d_perp_vec = [3 0 0].');

% All test molecules are planar in xy with normals equivalent up to sign.
assert(all(abs(T.normal_angle_deg) < 1e-10), ...
    'All molecule normals should be parallel in this descriptor test.');

descWithRef = builder.complete_molecule_descriptors_relative_to_reference( ...
    sys, 1, ...
    'StackAxis', [0 1 0], ...
    'IncludeReference', true, ...
    'IncludeNormals', false);

assert(ismember(1, descWithRef.table.molecule_id), ...
    'IncludeReference=true should include the reference molecule.');

rowRef = find(descWithRef.table.molecule_id == 1, 1, 'first');

assert(abs(descWithRef.table.distance(rowRef)) < 1e-12, ...
    'Reference descriptor row should have zero distance.');

end

function sys = local_make_descriptor_test_sys()

com = [
    5 5 5   % 1 reference
    5 7 5   % 2 same-stack +
    5 3 5   % 3 same-stack -
    8 5 5   % 4 side-stack
    2 5 5   % 5 side-stack opposite
    5 5 8   % 6 incomplete, should be excluded
];

isComplete = [true; true; true; true; true; false];

sys = local_make_planar_molecule_sys(com, isComplete);
sys.super_lattice = 10 * eye(3);

end

function sys = local_make_planar_molecule_sys(com, isComplete)

nMol = size(com, 1);
nPerMol = 4;
nSites = nMol * nPerMol;

local = [
    -0.5 -0.5 0
     0.5 -0.5 0
     0.5  0.5 0
    -0.5  0.5 0
];

site_pos = zeros(nSites, 3);
site_mol_id = zeros(nSites, 1);

site = 0;
site_indices = cell(nMol, 1);

for m = 1:nMol
    rows = site + (1:nPerMol);
    site_indices{m} = rows(:);

    site_pos(rows, :) = local + com(m, :);
    site_mol_id(rows) = m;

    site = site + nPerMol;
end

T = struct();
T.molecule_id = (1:nMol).';
T.unique_mol_id = T.molecule_id;
T.site_indices = site_indices;
T.n_sites = nPerMol * ones(nMol, 1);
T.com = com;
T.is_complete_in_display = isComplete(:);

sys = struct();
sys.site_pos = site_pos;
sys.site_mol_id = site_mol_id;
sys.unique_mol_id = site_mol_id;
sys.site_type = repmat({'C'}, nSites, 1);
sys.units.length = 'bohr';
sys.molecule_table = T;

end