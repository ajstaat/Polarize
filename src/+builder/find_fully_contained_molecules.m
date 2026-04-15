function Tinside = find_fully_contained_molecules(filename, sys, supercellSize, varargin)
%FIND_FULLY_CONTAINED_MOLECULES
% Find molecule images whose whole unwrapped geometry lies fully inside
% the supercell, using source CONTCAR molecules + image shifts from sys.
%
% Inputs
%   filename      : CONTCAR/POSCAR path
%   sys           : built Polarize system
%   supercellSize : [nx ny nz]
%
% Name-value options
%   'BondScale'   : default 1.20
%   'SortMolecules': default false
%   'Margin'      : default 0.05  (in supercell fractional coordinates)
%
% Output
%   Tinside       : subset of builder.list_molecules(sys) with fields:
%                   + minFrac, maxFrac, isFullyContained

    p = inputParser;
    addRequired(p, 'filename', @(x) ischar(x) || isstring(x));
    addRequired(p, 'sys', @isstruct);
    addRequired(p, 'supercellSize', @(x) isnumeric(x) && numel(x) == 3);
    addParameter(p, 'BondScale', 1.20, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'SortMolecules', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'Margin', 0.05, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    parse(p, filename, sys, supercellSize, varargin{:});
    opt = p.Results;

    % Unwrapped base molecules from source CONTCAR
    mols = io.unwrap_all_contcar_molecules(filename, ...
        'BondScale', opt.BondScale, ...
        'SortMolecules', opt.SortMolecules);

    T = builder.list_molecules(sys);
    n = height(T);

    minFrac = zeros(n, 3);
    maxFrac = zeros(n, 3);
    isInside = false(n, 1);

    dims = supercellSize(:).';   % [nx ny nz]

    for r = 1:n
        baseID = T.base_mol_id(r);
        shift = [T.ix(r), T.iy(r), T.iz(r)];

        Mol = mols{baseID};

        % Base molecule unwrapped fractional coords in unit-cell basis
        frac0 = Mol.fracUnwrapped;

        % Place into the supercell image
        fracImg = frac0 + shift;

        minFrac(r,:) = min(fracImg, [], 1);
        maxFrac(r,:) = max(fracImg, [], 1);

        isInside(r) = all(minFrac(r,:) >= opt.Margin) && ...
                      all(maxFrac(r,:) <= (dims - opt.Margin));
    end

    T.minFrac = minFrac;
    T.maxFrac = maxFrac;
    T.isFullyContained = isInside;

    Tinside = T(T.isFullyContained, :);
end