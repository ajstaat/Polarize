% src/+io/make_site_labels.m
function site_label = make_site_labels(species)
%MAKE_SITE_LABELS Generate labels like C1, C2, N1, H1...

    species = species(:);
    site_label = cell(size(species));
    counts = containers.Map('KeyType','char', 'ValueType','double');

    for i = 1:numel(species)
        sym = char(species{i});
        if isKey(counts, sym)
            counts(sym) = counts(sym) + 1;
        else
            counts(sym) = 1;
        end
        site_label{i} = sprintf('%s%d', sym, counts(sym));
    end
end