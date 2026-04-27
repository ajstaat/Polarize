function site_class = infer_simple_site_classes(species, A)
%INFER_SIMPLE_SITE_CLASSES Infer simple graph-based site classes.
%
% site_class = io.infer_simple_site_classes(species, A)
%
% Inputs
%   species   N x 1 cellstr/string array of element symbols
%   A         N x N adjacency matrix, logical/numeric/sparse
%
% Output
%   site_class N x 1 cell array of class labels:
%       C_deg3
%       C_deg4
%       H_on_C_deg3
%       H_on_C_deg4
%       N
%       O
%
% Notes
%   - This is intentionally simple and tailored to the current molecules.
%   - H is classified from the coordination of its bonded carbon neighbor.
%   - N and O are left as plain element classes.
%   - Unsupported local environments raise an error so they do not silently
%     get mistyped.

if isstring(species)
    species = cellstr(species(:));
elseif iscellstr(species)
    species = species(:);
else
    error('io:infer_simple_site_classes:BadSpecies', ...
        'species must be a string array or cell array of character vectors.');
end

n = numel(species);

if ~isequal(size(A), [n n])
    error('io:infer_simple_site_classes:BadAdjacency', ...
        'Adjacency matrix A must be N x N where N = numel(species).');
end

A = logical(A);
deg = full(sum(A, 2));

site_class = cell(n, 1);

for i = 1:n
    sym = strtrim(char(species{i}));
    nbrs = find(A(i, :));

    switch upper(sym)
        case 'C'
            switch deg(i)
                case 3
                    site_class{i} = 'C_deg3';

                case 4
                    site_class{i} = 'C_deg4';

                otherwise
                    error('io:infer_simple_site_classes:UnsupportedCarbonCoordination', ...
                        'Carbon at site %d has degree %d; expected 3 or 4.', ...
                        i, deg(i));
            end

        case 'H'
            if numel(nbrs) ~= 1
                error('io:infer_simple_site_classes:UnsupportedHydrogenCoordination', ...
                    'Hydrogen at site %d has %d bonded neighbors; expected 1.', ...
                    i, numel(nbrs));
            end

            j = nbrs(1);
            nbrSym = strtrim(char(species{j}));

            if ~strcmpi(nbrSym, 'C')
                error('io:infer_simple_site_classes:UnsupportedHydrogenNeighbor', ...
                    'Hydrogen at site %d is bonded to %s; only H-on-C is currently supported.', ...
                    i, nbrSym);
            end

            switch deg(j)
                case 3
                    site_class{i} = 'H_on_C_deg3';

                case 4
                    site_class{i} = 'H_on_C_deg4';

                otherwise
                    error('io:infer_simple_site_classes:UnsupportedHydrogenCarbonNeighbor', ...
                        ['Hydrogen at site %d is bonded to carbon of degree %d; ' ...
                         'expected 3 or 4.'], ...
                        i, deg(j));
            end

        case 'N'
            site_class{i} = 'N';

        case 'O'
            site_class{i} = 'O';

        otherwise
            error('io:infer_simple_site_classes:UnsupportedElement', ...
                'Unsupported element at site %d: %s', i, sym);
    end
end

end