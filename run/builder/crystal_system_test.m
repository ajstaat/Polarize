%% run_test_make_crystal_system
clear; clc;

rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
filename   = fullfile(rootFolder, 'a_0.0_CONTCAR.vasp');

%% ------------------------------------------------------------------------
% Import crystal template
crystal = io.import_contcar_as_crystal(filename, ...
    'BondScale', 1.20, ...
    'WrapFractional', true, ...
    'SortMolecules', false);

fprintf('Imported crystal template:\n');
fprintf('  nSites     = %d\n', size(crystal.cart_coords, 1));
fprintf('  nBaseMols  = %d\n', numel(unique(crystal.base_mol_id)));

%% ------------------------------------------------------------------------
% Minimal model
model = struct();
model.polarizable_types = {'H','C','N','O'};
model.alpha_by_type = struct( ...
    'H', 0.696, ...
    'C', 1.750, ...
    'N', 1.073, ...
    'O', 0.837);
model.thole_a = 0.39;

%% ------------------------------------------------------------------------
% Working options
opts = struct();
opts.supercell_size = [2 4 1];
opts.bondScale = 1.20;
opts.verbose = true;
opts.removeMolIDs = [];
opts.activeMolIDs = [];

%% ------------------------------------------------------------------------
% Build working system
sys = builder.make_crystal_system(crystal, model, opts);

fprintf('\nBuilt working system:\n');
fprintf('  n_unit_sites = %d\n', sys.n_unit_sites);
fprintf('  n_cells      = %d\n', sys.n_cells);
fprintf('  n_sites      = %d\n', sys.n_sites);

expectedSites = size(crystal.cart_coords,1) * prod(opts.supercell_size);
fprintf('  expected     = %d\n', expectedSites);

if sys.n_sites ~= expectedSites
    error('Site count mismatch: got %d, expected %d.', sys.n_sites, expectedSites);
end

%% ------------------------------------------------------------------------
% Basic required fields
requiredFields = { ...
    'site_pos', ...
    'site_frac', ...
    'site_type', ...
    'site_label', ...
    'base_mol_id', ...
    'site_mol_id', ...
    'molecule_table', ...
    'supercell_pbc_bond_graph'};

for k = 1:numel(requiredFields)
    name = requiredFields{k};
    if ~isfield(sys, name) || isempty(sys.(name))
        error('Missing or empty required field: sys.%s', name);
    end
end

fprintf('\nRequired fields present.\n');

%% ------------------------------------------------------------------------
% Site-level size checks
if size(sys.site_pos,1) ~= sys.n_sites || size(sys.site_pos,2) ~= 3
    error('sys.site_pos size mismatch.');
end

if size(sys.site_frac,1) ~= sys.n_sites || size(sys.site_frac,2) ~= 3
    error('sys.site_frac size mismatch.');
end

if numel(sys.site_type) ~= sys.n_sites
    error('sys.site_type size mismatch.');
end

if numel(sys.site_label) ~= sys.n_sites
    error('sys.site_label size mismatch.');
end

if numel(sys.base_mol_id) ~= sys.n_sites
    error('sys.base_mol_id size mismatch.');
end

if numel(sys.site_mol_id) ~= sys.n_sites
    error('sys.site_mol_id size mismatch.');
end

fprintf('Site-level array sizes are consistent.\n');

%% ------------------------------------------------------------------------
% Molecule-table checks
T = sys.molecule_table;

tableFields = { ...
    'molecule_id', ...
    'site_indices', ...
    'n_sites', ...
    'com', ...
    'is_complete_in_display', ...
    'n_display_fragments', ...
    'largest_fragment_fraction'};

for k = 1:numel(tableFields)
    name = tableFields{k};
    if ~isfield(T, name) || isempty(T.(name))
        error('sys.molecule_table.%s is missing or empty.', name);
    end
end

nMol = numel(T.molecule_id);

fprintf('\nSupercell molecule table:\n');
fprintf('  n_molecules = %d\n', nMol);

if numel(T.site_indices) ~= nMol
    error('molecule_table.site_indices size mismatch.');
end

if numel(T.n_sites) ~= nMol
    error('molecule_table.n_sites size mismatch.');
end

if size(T.com,1) ~= nMol || size(T.com,2) ~= 3
    error('molecule_table.com size mismatch.');
end

if numel(T.is_complete_in_display) ~= nMol
    error('molecule_table.is_complete_in_display size mismatch.');
end

if numel(T.n_display_fragments) ~= nMol
    error('molecule_table.n_display_fragments size mismatch.');
end

if numel(T.largest_fragment_fraction) ~= nMol
    error('molecule_table.largest_fragment_fraction size mismatch.');
end

%% ------------------------------------------------------------------------
% Molecule-table consistency with site_mol_id
allMolIDs = unique(sys.site_mol_id(:), 'stable');

if numel(allMolIDs) ~= nMol
    error('Mismatch between unique(sys.site_mol_id) and molecule_table rows.');
end

if ~isequal(allMolIDs(:), T.molecule_id(:))
    error('molecule_table.molecule_id does not match unique(sys.site_mol_id).');
end

for m = 1:nMol
    idx = T.site_indices{m};
    molID = T.molecule_id(m);

    if isempty(idx)
        error('Molecule %d has empty site_indices.', molID);
    end

    if ~all(sys.site_mol_id(idx) == molID)
        error('site_indices mismatch for molecule %d.', molID);
    end

    if numel(idx) ~= T.n_sites(m)
        error('n_sites mismatch for molecule %d.', molID);
    end
end

fprintf('Molecule-table site membership is consistent.\n');

%% ------------------------------------------------------------------------
% Completeness summary
nComplete = nnz(T.is_complete_in_display);
fprintf('\nDisplayed-supercell completeness:\n');
fprintf('  complete molecules = %d / %d\n', nComplete, nMol);

disp(table( ...
    T.molecule_id(:), ...
    T.n_sites(:), ...
    T.is_complete_in_display(:), ...
    T.n_display_fragments(:), ...
    T.largest_fragment_fraction(:), ...
    'VariableNames', { ...
        'molecule_id', ...
        'n_sites', ...
        'is_complete', ...
        'n_fragments', ...
        'largest_fragment_fraction'}))

if nComplete == 0
    warning('No complete molecules found in displayed supercell.');
end

%% ------------------------------------------------------------------------
% Pick centered reference from complete molecules
completeRows = find(T.is_complete_in_display);

if isempty(completeRows)
    warning('No complete molecules available for centered-reference selection.');
    refMolID = [];
else
    center = 0.5 * sum(sys.super_lattice, 1);

    comComplete = T.com(completeRows, :);
    d2 = sum((comComplete - center).^2, 2);
    [~, iMin] = min(d2);

    chosenRow = completeRows(iMin);
    refMolID = T.molecule_id(chosenRow);

    fprintf('\nCentered complete reference molecule:\n');
    fprintf('  molecule_id = %d\n', refMolID);
    fprintf('  COM         = [%9.4f %9.4f %9.4f]\n', T.com(chosenRow, :));
    fprintf('  center      = [%9.4f %9.4f %9.4f]\n', center);
end

%% ------------------------------------------------------------------------
% Active molecule smoke test
if ~isempty(refMolID)
    sys2 = builder.select_active_molecules(sys, refMolID);

    fprintf('\nActive molecule test:\n');
    fprintf('  requested active ID = %d\n', refMolID);
    fprintf('  active sites        = %d\n', nnz(sys2.site_is_active));

    idx = builder.site_indices_for_molecule(sys2, refMolID);
    if ~all(sys2.site_is_active(idx))
        error('Active mask failed for molecule %d.', refMolID);
    end
end

%% ------------------------------------------------------------------------
% Remove molecule smoke test
if nMol >= 2
    removeID = T.molecule_id(2);
    sys3 = builder.remove_molecules(sys, removeID);

    fprintf('\nRemove molecule test:\n');
    fprintf('  removed ID  = %d\n', removeID);
    fprintf('  new n_sites = %d\n', sys3.n_sites);

    if any(sys3.site_mol_id == removeID)
        error('Removed molecule ID still present in site_mol_id.');
    end
end

fprintf('\nCheckpoint passed.\n');