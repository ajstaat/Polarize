function molecules = unwrap_all_contcar_molecules(inputStructOrFilename, varargin)
%UNWRAP_ALL_CONTCAR_MOLECULES
% Identify molecules in a periodic unit cell and unwrap each molecule into
% a coherent Cartesian representation.
%
% Input
%   inputStructOrFilename : either
%       - filename for a CONTCAR/POSCAR-style file, or
%       - struct S with fields:
%             .species
%             .frac
%             .cart   (optional; recomputed from frac/lattice if absent)
%             .lattice
%
% Optional name-value inputs
%   'BondScale'     : scalar, default 1.20
%   'SortMolecules' : logical, default false
%
% Output
%   molecules{k} is a struct with fields:
%       .indices
%       .labels
%       .fracWrapped
%       .fracUnwrapped
%       .cart
%       .com
%       .base_mol_id
%
% Notes
%   - Molecules are identified from a PBC-aware bond graph on the UNIT CELL.
%   - This is intended to build the unit-cell chemical template, not final
%     molecule identities for a built supercell.

    p = inputParser;
    addRequired(p, 'inputStructOrFilename', @(x) ...
        (ischar(x) || isstring(x)) || isstruct(x));
    addParameter(p, 'BondScale', 1.20, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'SortMolecules', false, @(x) islogical(x) && isscalar(x));
    parse(p, inputStructOrFilename, varargin{:});
    opt = p.Results;

    % -------------------------------------------------------------
    % Accept either a filename or a pre-parsed structure
    % -------------------------------------------------------------
    if isstruct(inputStructOrFilename)
        S = inputStructOrFilename;
        requiredFields = {'species','frac','lattice'};
        for k = 1:numel(requiredFields)
            if ~isfield(S, requiredFields{k})
                error('io:unwrap_all_contcar_molecules:BadStructureInput', ...
                    'Input struct is missing required field: %s', requiredFields{k});
            end
        end

        if ~isfield(S, 'cart') || isempty(S.cart)
            S.cart = S.frac * S.lattice;
        end
    else
        S = io.read_vasp_structure(inputStructOrFilename);
    end

    % -------------------------------------------------------------
    % Normalize wrapped unit-cell coordinates
    % -------------------------------------------------------------
    S.frac = S.frac - floor(S.frac);
    S.cart = S.frac * S.lattice;

    % -------------------------------------------------------------
    % PBC-aware molecule identification on the UNIT CELL
    % -------------------------------------------------------------
    [A, ~] = io.build_pbc_bond_graph(S.species, S.frac, S.lattice, opt.BondScale);

    componentID = conncomp(graph(A)).';
    molIDs = unique(componentID, 'stable');

    molecules = cell(numel(molIDs), 1);
    for k = 1:numel(molIDs)
        idx = find(componentID == molIDs(k));
        mol = io.unwrap_vasp_molecule(S, idx, A);
        mol.base_mol_id = k;
        molecules{k} = mol;
    end

    % -------------------------------------------------------------
    % Optional sorting by molecule COM
    % -------------------------------------------------------------
    if opt.SortMolecules
        coms = cellfun(@(m) m.com, molecules, 'UniformOutput', false);
        comMat = vertcat(coms{:});
        [~, order] = sortrows(comMat, [3 2 1]);
        molecules = molecules(order);

        % Renumber after sorting so base_mol_id matches output order
        for k = 1:numel(molecules)
            molecules{k}.base_mol_id = k;
        end
    end
end