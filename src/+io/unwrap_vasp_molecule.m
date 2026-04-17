function Mol = unwrap_vasp_molecule(S, atomIdx, A)
%UNWRAP_VASP_MOLECULE Unwrap one molecular fragment under PBC.
%
% Mol = io.unwrap_vasp_molecule(S, atomIdx)
% Mol = io.unwrap_vasp_molecule(S, atomIdx, A)
%
% Inputs
%   S       struct with fields:
%             .species
%             .frac
%             .lattice
%           and optionally:
%             .cart
%
%   atomIdx  vector of site indices belonging to one connected component
%
%   A        optional full adjacency matrix for the structure
%
% Output
%   Mol struct with fields:
%       .indices
%       .labels
%       .fracWrapped
%       .fracUnwrapped
%       .cart
%       .com
%
% Notes
%   - The input atomIdx is assumed to define one molecule in the periodic
%     unit cell, typically from a PBC-aware bond graph.
%   - Unwrapping is done by traversing the bonded graph and applying
%     minimum-image fractional displacements.
%   - The order of atoms is preserved from atomIdx.

    atomIdx = atomIdx(:);
    nloc = numel(atomIdx);

    if nloc == 0
        error('io:unwrap_vasp_molecule:EmptyAtomIdx', ...
            'atomIdx must contain at least one site index.');
    end

    validate_structure_input(S);

    if any(atomIdx < 1) || any(atomIdx > numel(S.species)) || any(mod(atomIdx,1) ~= 0)
        error('io:unwrap_vasp_molecule:BadAtomIdx', ...
            'atomIdx must contain valid integer site indices.');
    end

    if nargin < 3 || isempty(A)
        [A, ~] = io.build_pbc_bond_graph(S.species, S.frac, S.lattice, 1.20);
    end

    if ~isequal(size(A), [numel(S.species), numel(S.species)])
        error('io:unwrap_vasp_molecule:BadAdjacency', ...
            'Adjacency matrix A must be N x N where N = numel(S.species).');
    end

    Aloc = A(atomIdx, atomIdx);

    fracWrapped = S.frac(atomIdx, :);
    labels = S.species(atomIdx);

    fracUnwrapped = nan(nloc, 3);
    visited = false(nloc, 1);

    % Anchor the first atom at its wrapped fractional position
    fracUnwrapped(1, :) = fracWrapped(1, :);
    visited(1) = true;
    queue = 1;

    % Breadth-first traversal over the local bonded graph
    while ~isempty(queue)
        current = queue(1);
        queue(1) = [];

        nbrs = find(Aloc(current, :));
        for nb = nbrs
            if ~visited(nb)
                deltaWrapped = fracWrapped(nb, :) - fracWrapped(current, :);
                deltaMin = deltaWrapped - round(deltaWrapped);

                fracUnwrapped(nb, :) = fracUnwrapped(current, :) + deltaMin;
                visited(nb) = true;
                queue(end+1) = nb; %#ok<AGROW>
            end
        end
    end

    % Fallback for disconnected leftovers (should not happen if atomIdx
    % came from a connected component, but keeps behavior robust)
    for i = 1:nloc
        if ~visited(i)
            fracUnwrapped(i, :) = fracWrapped(i, :);
        end
    end

    Mol = struct();
    Mol.indices = atomIdx;
    Mol.labels = labels;
    Mol.fracWrapped = fracWrapped;
    Mol.fracUnwrapped = fracUnwrapped;
    Mol.cart = fracUnwrapped * S.lattice;
    Mol.com = mean(Mol.cart, 1);
end


function validate_structure_input(S)
    requiredFields = {'species', 'frac', 'lattice'};
    for k = 1:numel(requiredFields)
        name = requiredFields{k};
        if ~isfield(S, name) || isempty(S.(name))
            error('io:unwrap_vasp_molecule:BadStructureInput', ...
                'Input struct S is missing required field: %s', name);
        end
    end

    species = S.species(:);
    n = numel(species);

    if size(S.frac, 1) ~= n || size(S.frac, 2) ~= 3
        error('io:unwrap_vasp_molecule:BadFracCoords', ...
            'S.frac must be N x 3 with N = numel(S.species).');
    end

    if ~isequal(size(S.lattice), [3 3])
        error('io:unwrap_vasp_molecule:BadLattice', ...
            'S.lattice must be 3 x 3.');
    end
end