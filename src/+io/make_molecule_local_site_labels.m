function site_label = make_molecule_local_site_labels(species, mol_id)
%MAKE_MOLECULE_LOCAL_SITE_LABELS Generate per-molecule local atom labels.
%
% site_label = io.make_molecule_local_site_labels(species, mol_id)
%
% Inputs
%   species   N x 1 cellstr/string array of element symbols
%   mol_id    N x 1 molecule identifiers
%
% Output
%   site_label N x 1 cell array of labels like C1, C2, H1, N1, ...
%
% Labels are counted separately within each molecule. For example, two
% molecules each containing carbon atoms will each start from C1.

if isstring(species)
    species = cellstr(species(:));
elseif iscellstr(species)
    species = species(:);
else
    error('io:make_molecule_local_site_labels:BadSpecies', ...
        'species must be a string array or cell array of character vectors.');
end

mol_id = mol_id(:);

if numel(species) ~= numel(mol_id)
    error('io:make_molecule_local_site_labels:SizeMismatch', ...
        'species and mol_id must have the same length.');
end

n = numel(species);
site_label = cell(n, 1);

mols = unique(mol_id, 'stable');

for m = 1:numel(mols)
    idx = find(mol_id == mols(m));

    counts = containers.Map('KeyType', 'char', 'ValueType', 'double');

    for k = 1:numel(idx)
        i = idx(k);

        sym = strtrim(char(species{i}));

        if isKey(counts, sym)
            counts(sym) = counts(sym) + 1;
        else
            counts(sym) = 1;
        end

        site_label{i} = sprintf('%s%d', sym, counts(sym));
    end
end

end