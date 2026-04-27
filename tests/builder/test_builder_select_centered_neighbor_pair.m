function test_builder_select_centered_neighbor_pair()
%TEST_BUILDER_SELECT_CENTERED_NEIGHBOR_PAIR Verify centered pair selection.

sys = local_make_centered_pair_sys();

result = builder.select_centered_neighbor_pair(sys, ...
    'Relation', 'same_stack', ...
    'Shell', 1, ...
    'StackAxis', [0 1 0], ...
    'Direction', 'either', ...
    'PerpTol', 1e-6, ...
    'ShellTol', 1e-6, ...
    'Verbose', false);

assert(strcmp(result.relation, 'same_stack'), ...
    'Unexpected relation label.');

assert(ismember(result.reference_mol_id, [1 2]), ...
    'Chosen reference should be one of the centered pair molecules.');

assert(ismember(result.neighbor_mol_id, [1 2]), ...
    'Chosen neighbor should be one of the centered pair molecules.');

assert(result.reference_mol_id ~= result.neighbor_mol_id, ...
    'Reference and neighbor molecule IDs should differ.');

assert(norm(result.pair_midpoint - [5 5 5]) < 1e-12, ...
    'Chosen pair midpoint should be exactly at the supercell center.');

assert(abs(result.midpoint_distance) < 1e-12, ...
    'Chosen pair midpoint distance should be zero.');

assert(~isempty(result.candidate_table), ...
    'Centered pair selector should return a candidate table.');

end

function sys = local_make_centered_pair_sys()

% Molecules 1 and 2 are same-stack around the center:
%   COM1 = [5 3 5]
%   COM2 = [5 7 5]
% midpoint = [5 5 5]
%
% Molecules 3 and 4 form another same-stack pair but off-center.

com = [
    5 3 5
    5 7 5
    2 3 5
    2 7 5
];

isComplete = true(size(com, 1), 1);

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