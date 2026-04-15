function molecules = unwrap_all_contcar_molecules(filename, varargin)
%UNWRAP_ALL_CONTCAR_MOLECULES
% Read a CONTCAR/POSCAR, identify molecules under PBC, and unwrap each
% molecule into a coherent Cartesian representation.
%
% Output:
%   molecules{k} is a struct with fields:
%       .indices
%       .labels
%       .fracWrapped
%       .fracUnwrapped
%       .cart
%       .com
%
% Example:
%   mols = unwrap_all_contcar_molecules('CONTCAR.vasp', ...
%       'BondScale', 1.20, 'SortMolecules', true);

    p = inputParser;
    addRequired(p, 'filename', @(x) ischar(x) || isstring(x));
    addParameter(p, 'BondScale', 1.20, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'SortMolecules', false, @(x) islogical(x) && isscalar(x));
    parse(p, filename, varargin{:});
    opt = p.Results;

    S = io.read_vasp_structure(filename);

    % Keep wrapped coordinates in [0,1)
    S.frac = S.frac - floor(S.frac);
    S.cart = S.frac * S.lattice;

    [A, ~] = io.build_vasp_bond_graph(S, opt.BondScale);
    comp = conncomp(graph(A)).';
    molIDs = unique(comp, 'stable');

    molecules = cell(numel(molIDs), 1);
    for k = 1:numel(molIDs)
        idx = find(comp == molIDs(k));
        molecules{k} = io.unwrap_vasp_molecule(S, idx);
    end

    if opt.SortMolecules
        coms = cellfun(@(m) m.com, molecules, 'UniformOutput', false);
        comMat = vertcat(coms{:});
        [~, order] = sortrows(comMat, [3 2 1]);
        molecules = molecules(order);
    end
end