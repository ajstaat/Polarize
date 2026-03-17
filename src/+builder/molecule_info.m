function info = molecule_info(sys, uniqueMolID)
%MOLECULE_INFO Return bookkeeping info for one unique molecule.

row = builder.find_molecule_row(sys, uniqueMolID);

info = struct();
info.unique_mol_id = sys.molecule_table.unique_mol_id(row);
info.base_mol_id = sys.molecule_table.base_mol_id(row);
info.cell_shift = sys.molecule_table.cell_shift(row, :);
info.site_indices = sys.molecule_table.site_indices{row};

if isfield(sys, 'site_pos') && ~isempty(sys.site_pos)
    coords = sys.site_pos(info.site_indices, :);
    info.center = mean(coords, 1);
else
    info.center = [];
end

end