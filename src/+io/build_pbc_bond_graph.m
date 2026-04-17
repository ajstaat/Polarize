function [A, Dcart] = build_pbc_bond_graph(species, frac_coords, lattice, scaleFactor, varargin)
%BUILD_PBC_BOND_GRAPH PBC-aware bond graph using minimum-image distances.
%
% Inputs
%   species      N x 1 cell array of element symbols
%   frac_coords  N x 3 fractional coordinates
%   lattice      3 x 3 lattice matrix, rows are lattice vectors
%   scaleFactor  scalar multiplier on covalent-radius sum, default 1.20
%
% Optional name-value inputs
%   'Verbose'           logical, default false
%   'ProgressInterval'  positive integer, default 250
%
% Outputs
%   A            N x N logical adjacency matrix
%   Dcart        N x N minimum-image Cartesian distance matrix

    if nargin < 4 || isempty(scaleFactor)
        scaleFactor = 1.20;
    end

    p = inputParser;
    addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ProgressInterval', 250, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    parse(p, varargin{:});
    opt = p.Results;

    species = species(:);
    n = numel(species);

    if size(frac_coords, 1) ~= n || size(frac_coords, 2) ~= 3
        error('io:build_pbc_bond_graph:BadFracCoords', ...
            'frac_coords must be N x 3 and match species.');
    end

    if ~isequal(size(lattice), [3 3])
        error('io:build_pbc_bond_graph:BadLattice', ...
            'lattice must be 3 x 3.');
    end

    if ~isscalar(scaleFactor) || ~isfinite(scaleFactor) || scaleFactor <= 0
        error('io:build_pbc_bond_graph:BadScaleFactor', ...
            'scaleFactor must be a positive finite scalar.');
    end

    if opt.Verbose
        fprintf('build_pbc_bond_graph: starting for %d atoms\n', n);
        tStart = tic;
    end

    A = false(n, n);
    Dcart = zeros(n, n);

    for i = 1:n-1
        si = species{i};

        for j = i+1:n
            sj = species{j};

            df = frac_coords(i,:) - frac_coords(j,:);
            df = df - round(df);
            dcart = norm(df * lattice);

            Dcart(i,j) = dcart;
            Dcart(j,i) = dcart;

            if io.bond_graph_tools.is_bonded(si, sj, dcart, scaleFactor)
                A(i,j) = true;
                A(j,i) = true;
            end
        end

        if opt.Verbose
            if mod(i, opt.ProgressInterval) == 0 || i == n-1
                fprintf('  processed %d / %d outer-loop atoms\n', i, n-1);
            end
        end
    end

    if opt.Verbose
        fprintf('build_pbc_bond_graph: finished in %.3f s\n', toc(tStart));
    end
end