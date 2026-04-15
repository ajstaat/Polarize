% src/+io/unwrap_vasp_molecule.m
function Mol = unwrap_vasp_molecule(S, atomIdx)
%UNWRAP_VASP_MOLECULE Unwrap one molecular fragment under PBC.
% atomIdx is a vector of site indices belonging to one connected component.

    atomIdx = atomIdx(:);
    nloc = numel(atomIdx);

    [A, ~] = io.build_vasp_bond_graph(S, 1.20);
    Aloc = A(atomIdx, atomIdx);

    fracWrapped = S.frac(atomIdx,:);
    labels = S.species(atomIdx);

    fracUnwrapped = nan(nloc,3);
    visited = false(nloc,1);

    fracUnwrapped(1,:) = fracWrapped(1,:);
    visited(1) = true;
    queue = 1;

    while ~isempty(queue)
        current = queue(1);
        queue(1) = [];

        nbrs = find(Aloc(current,:));
        for nb = nbrs
            if ~visited(nb)
                deltaWrapped = fracWrapped(nb,:) - fracWrapped(current,:);
                deltaMin = deltaWrapped - round(deltaWrapped);
                fracUnwrapped(nb,:) = fracUnwrapped(current,:) + deltaMin;
                visited(nb) = true;
                queue(end+1) = nb; %#ok<AGROW>
            end
        end
    end

    for i = 1:nloc
        if ~visited(i)
            fracUnwrapped(i,:) = fracWrapped(i,:);
        end
    end

    Mol = struct();
    Mol.indices = atomIdx;
    Mol.labels  = labels;
    Mol.fracWrapped   = fracWrapped;
    Mol.fracUnwrapped = fracUnwrapped;
    Mol.cart = fracUnwrapped * S.lattice;
    Mol.com  = mean(Mol.cart, 1);
end