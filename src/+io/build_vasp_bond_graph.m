% src/+io/build_vasp_bond_graph.m
function [A, Dcart] = build_vasp_bond_graph(S, scaleFactor)
%BUILD_VASP_BOND_GRAPH PBC-aware bond graph using minimum-image distances.

    n = numel(S.species);
    A = false(n,n);
    Dcart = zeros(n,n);

    for i = 1:n-1
        ri = covalent_radius(S.species{i});
        for j = i+1:n
            rj = covalent_radius(S.species{j});

            df = S.frac(i,:) - S.frac(j,:);
            df = df - round(df);
            dcart = norm(df * S.lattice);

            Dcart(i,j) = dcart;
            Dcart(j,i) = dcart;

            cutoff = scaleFactor * (ri + rj);

            if strcmpi(S.species{i}, 'H') && strcmpi(S.species{j}, 'H')
                continue;
            end

            if dcart < cutoff
                A(i,j) = true;
                A(j,i) = true;
            end
        end
    end
end

function r = covalent_radius(sym)
    switch upper(sym)
        case 'H'
            r = 0.31;
        case 'C'
            r = 0.76;
        case 'N'
            r = 0.71;
        case 'O'
            r = 0.66;
        case 'F'
            r = 0.57;
        case 'S'
            r = 1.05;
        case 'CL'
            r = 1.02;
        case 'BR'
            r = 1.20;
        case 'I'
            r = 1.39;
        case 'SI'
            r = 1.11;
        otherwise
            error('No covalent radius defined for element: %s', sym);
    end
end