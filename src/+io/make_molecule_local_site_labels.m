function site_label = make_molecule_local_site_labels(species, mol_id)
%MAKE_MOLECULE_LOCAL_SITE_LABELS
% Generate labels like C1, C2, N1... separately within each molecule.

    species = species(:);
    mol_id = mol_id(:);

    if numel(species) ~= numel(mol_id)
        error('species and mol_id must have the same length.');
    end

    site_label = cell(size(species));
    mols = unique(mol_id, 'stable');

    for m = 1:numel(mols)
        idx = find(mol_id == mols(m));

        counts = containers.Map('KeyType','char', 'ValueType','double');

        for k = 1:numel(idx)
            i = idx(k);
            sym = char(species{i});

            if isKey(counts, sym)
                counts(sym) = counts(sym) + 1;
            else
                counts(sym) = 1;
            end

            site_label{i} = sprintf('%s%d', sym, counts(sym));
        end
    end
end