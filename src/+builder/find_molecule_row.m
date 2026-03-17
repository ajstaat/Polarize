function row = find_molecule_row(sys, uniqueMolID)
%FIND_MOLECULE_ROW Find row index in molecule_table for a given unique molecule ID.

if ~isfield(sys, 'molecule_table') || isempty(sys.molecule_table)
    error('sys.molecule_table is missing or empty.');
end

row = find(sys.molecule_table.unique_mol_id == uniqueMolID);

if isempty(row)
    error('uniqueMolID %d not found in molecule_table.', uniqueMolID);
end

if numel(row) > 1
    error('uniqueMolID %d appears multiple times in molecule_table.', uniqueMolID);
end

end