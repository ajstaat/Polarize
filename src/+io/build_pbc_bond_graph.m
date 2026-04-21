function [A, Dcart] = build_pbc_bond_graph(species, frac_coords, lattice, scaleFactor, varargin)
%BUILD_PBC_BOND_GRAPH PBC-aware bond graph using shared cell-list search.
%
% [A, Dcart] = io.build_pbc_bond_graph(species, frac_coords, lattice, scaleFactor)
% [A, Dcart] = io.build_pbc_bond_graph(..., 'Verbose', tf, ...
%                                           'ProgressInterval', k, ...
%                                           'ReturnDcart', tf)
%
% Inputs
%   species      N x 1 cell array of element symbols
%   frac_coords  N x 3 fractional coordinates
%   lattice      3 x 3 lattice matrix, ROWS are lattice vectors
%   scaleFactor  scalar multiplier on covalent-radius sum, default 1.20
%
% Optional name-value inputs
%   'Verbose'          logical, default false
%   'ProgressInterval' positive integer, default 2500
%   'ReturnDcart'      logical, default false
%
% Outputs
%   A      N x N logical adjacency matrix
%   Dcart  N x N minimum-image Cartesian distance matrix if requested;
%          otherwise []
%
% Notes
% - Candidate bonded pairs are generated through the shared
%   geom.build_spatial_index / geom.query_pairs_within_cutoff interface.
% - Exact species-specific bond tests remain delegated to
%   io.bond_graph_tools.is_bonded(...).
% - lattice is assumed to have lattice vectors as ROWS, consistent with
%   existing io conventions.
% - geom periodic routines expect lattice vectors as COLUMNS, so we pass
%   lattice.' into geom.build_spatial_index(...).

if nargin < 4 || isempty(scaleFactor)
    scaleFactor = 1.20;
end

p = inputParser;
addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ProgressInterval', 2500, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'ReturnDcart', false, @(x) islogical(x) && isscalar(x));
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

if opt.ReturnDcart
    Dcart = zeros(n, n);
else
    Dcart = [];
end

% Precompute covalent radii once.
radii = zeros(n, 1);
for i = 1:n
    radii(i) = io.bond_graph_tools.covalent_radius(species{i});
end

% Convert fractional coordinates to Cartesian using row-vector convention.
cart_coords = frac_coords * lattice;

% Global candidate cutoff: any bond distance must be <= scaleFactor *
% (r_cov(i) + r_cov(j)), so a safe global bound is scaleFactor * 2*max(radii).
maxBondCutoff = scaleFactor * (2 * max(radii));

% geom periodic routines expect lattice vectors as columns.
cellCols = lattice.';

spatial = geom.build_spatial_index(cart_coords, struct( ...
    'isPeriodic', true, ...
    'cell', cellCols, ...
    'method', 'cell_list', ...
    'cutoff', maxBondCutoff));

pairs = geom.query_pairs_within_cutoff(spatial, maxBondCutoff, struct( ...
    'return_r', true, ...
    'return_dr', false, ...
    'return_full_idx', true));

nPairs = numel(pairs.i);

if opt.Verbose
    fprintf('build_pbc_bond_graph: candidate pairs within global cutoff = %d\n', nPairs);
end

for kk = 1:nPairs
    i = pairs.i(kk);
    j = pairs.j(kk);
    dij = pairs.r(kk);

    if io.bond_graph_tools.is_bonded(species{i}, species{j}, dij, scaleFactor)
        A(i, j) = true;
        A(j, i) = true;
    end

    if opt.Verbose
        if mod(kk, opt.ProgressInterval) == 0 || kk == nPairs
            fprintf(' processed %d / %d candidate bond pairs\n', kk, nPairs);
        end
    end
end

if opt.ReturnDcart
    % Dense all-pairs minimum-image distance matrix, only when explicitly requested.
    for i = 1:(n-1)
        df = frac_coords((i+1):n, :) - frac_coords(i, :);
        df = df - round(df);
        dcart = df * lattice;
        d = sqrt(sum(dcart.^2, 2));
        Dcart(i, (i+1):n) = d.';
        Dcart((i+1):n, i) = d.';
    end
end

if opt.Verbose
    fprintf('build_pbc_bond_graph: finished in %.3f s\n', toc(tStart));
end

end