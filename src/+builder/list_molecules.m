function T = list_molecules(sys)
%LIST_MOLECULES Build a summary table of molecules in the working system.
%
% Output
%   T   MATLAB table with one row per unique molecule

if ~isfield(sys, 'molecule_table') || isempty(sys.molecule_table)
    error('sys.molecule_table is missing or empty.');
end

nMol = numel(sys.molecule_table.unique_mol_id);

unique_mol_id = sys.molecule_table.unique_mol_id(:);
base_mol_id   = sys.molecule_table.base_mol_id(:);

ix = sys.molecule_table.cell_shift(:, 1);
iy = sys.molecule_table.cell_shift(:, 2);
iz = sys.molecule_table.cell_shift(:, 3);

n_sites = zeros(nMol, 1);
cx = nan(nMol, 1);
cy = nan(nMol, 1);
cz = nan(nMol, 1);

for m = 1:nMol
    idx = sys.molecule_table.site_indices{m};
    n_sites(m) = numel(idx);

    if isfield(sys, 'site_pos') && ~isempty(sys.site_pos)
        ctr = mean(sys.site_pos(idx, :), 1);
        cx(m) = ctr(1);
        cy(m) = ctr(2);
        cz(m) = ctr(3);
    end
end

T = table(unique_mol_id, base_mol_id, ix, iy, iz, n_sites, cx, cy, cz);

end