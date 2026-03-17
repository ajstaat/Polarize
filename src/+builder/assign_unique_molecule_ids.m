function [unique_mol_id, mol_table] = assign_unique_molecule_ids(base_mol_id, cell_shift)
%ASSIGN_UNIQUE_MOLECULE_IDS Assign one unique molecule ID to each molecule image.
%
% Inputs
%   base_mol_id : N x 1 vector
%                 molecule ID in the unit cell for each site
%
%   cell_shift  : N x 3 integer translation vector [ix iy iz] for each site
%
% Outputs
%   unique_mol_id : N x 1 vector
%                   unique molecule ID for each site in the supercell
%
%   mol_table     : struct with one row per unique molecule:
%                   .unique_mol_id
%                   .base_mol_id
%                   .cell_shift
%                   .site_indices   (cell array)

if size(base_mol_id, 2) ~= 1
    base_mol_id = base_mol_id(:);
end

nSites = numel(base_mol_id);

if size(cell_shift, 1) ~= nSites || size(cell_shift, 2) ~= 3
    error('cell_shift must be N x 3 and match base_mol_id.');
end

% Key per site: [base_mol_id, ix, iy, iz]
keys = [base_mol_id, cell_shift];

% unique rows tells us which sites belong to the same molecule image
[uniqueKeys, ~, ic] = unique(keys, 'rows', 'stable');

nMol = size(uniqueKeys, 1);
unique_mol_id = ic;

% Build molecule table
mol_table = struct();
mol_table.unique_mol_id = (1:nMol).';
mol_table.base_mol_id = uniqueKeys(:, 1);
mol_table.cell_shift = uniqueKeys(:, 2:4);
mol_table.site_indices = cell(nMol, 1);

for m = 1:nMol
    mol_table.site_indices{m} = find(ic == m);
end

end