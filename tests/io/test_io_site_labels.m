function test_io_site_labels()
%TEST_IO_SITE_LABELS Verify per-molecule local atom labels.
%
% Labels should count separately inside each molecule.

species_one = {
    'C'
    'C'
    'C'
    'H'
    'H'
    'H'
    'H'
    'H'
    'H'
};

species = [species_one; species_one];

mol_id = [
    ones(numel(species_one), 1)
    2 * ones(numel(species_one), 1)
];

labels = io.make_molecule_local_site_labels(species, mol_id);

expected_one = {
    'C1'
    'C2'
    'C3'
    'H1'
    'H2'
    'H3'
    'H4'
    'H5'
    'H6'
};

expected = [expected_one; expected_one];

assert(isequal(labels, expected), ...
    'Local site labels should restart for each molecule.');

% Also check string species input.
labels_string = io.make_molecule_local_site_labels(string(species), mol_id);

assert(isequal(labels_string, expected), ...
    'String species input should produce the same labels.');

end