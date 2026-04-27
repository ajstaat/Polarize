function test_builder_choose_center_reference_molecule()
%TEST_BUILDER_CHOOSE_CENTER_REFERENCE_MOLECULE Verify center-reference choice.

sys = local_make_center_test_sys();

[refMolID, summary] = builder.choose_center_reference_molecule(sys, ...
    'RequireComplete', true, ...
    'Verbose', false);

assert(refMolID == 2, ...
    'Expected molecule 2 to be the complete molecule closest to the supercell center.');

assert(isequal(summary.chosen_mol_id, refMolID), ...
    'summary.chosen_mol_id should match returned refMolID.');

assert(isequal(summary.supercell_center, [5 5 5]), ...
    'Unexpected supercell center.');

assert(~ismember(3, summary.candidate_mol_ids), ...
    'Incomplete molecule 3 should not be considered when RequireComplete=true.');

[refAny, summaryAny] = builder.choose_center_reference_molecule(sys, ...
    'RequireComplete', false, ...
    'Verbose', false);

assert(refAny == 3, ...
    'With RequireComplete=false, molecule 3 at the exact center should be chosen.');

assert(ismember(3, summaryAny.candidate_mol_ids), ...
    'Incomplete molecule 3 should be considered when RequireComplete=false.');

end

function sys = local_make_center_test_sys()

sys = struct();
sys.super_lattice = 10 * eye(3);

T = struct();
T.molecule_id = [1; 2; 3];
T.unique_mol_id = T.molecule_id;
T.com = [
    1 1 1
    4 5 5
    5 5 5
];
T.n_sites = [4; 4; 4];
T.is_complete_in_display = [true; true; false];

sys.molecule_table = T;

end