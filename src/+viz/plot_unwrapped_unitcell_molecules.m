function out = plot_unwrapped_unitcell_molecules(input, varargin)
%PLOT_UNWRAPPED_UNITCELL_MOLECULES Plot molecules unwrapped from a POSCAR/CONTCAR.
%
% out = viz.plot_unwrapped_unitcell_molecules(filename)
% out = viz.plot_unwrapped_unitcell_molecules(S)
%
% Options
%   'Axes'           axes handle
%   'BondScale'      default 1.20
%   'SortMolecules'  default false
%   'DrawBox'        default true
%   'ShowBonds'      default true
%   'MarkerSize'     default 60
%   'Title'          default ''

p = inputParser;
addRequired(p, 'input');
addParameter(p, 'Axes', [], @(x) isempty(x) || isgraphics(x, 'axes'));
addParameter(p, 'BondScale', 1.20, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'SortMolecules', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'DrawBox', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ShowBonds', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'MarkerSize', 60, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));
parse(p, input, varargin{:});

if ischar(input) || isstring(input)
    S = io.read_vasp_structure(input);
elseif isstruct(input)
    S = input;
else
    error('viz:plot_unwrapped_unitcell_molecules:BadInput', ...
        'Input must be a filename or parsed VASP struct.');
end

[A, ~] = io.build_pbc_bond_graph(S.species, S.frac, S.lattice, p.Results.BondScale);

molecules = io.unwrap_all_contcar_molecules(S, ...
    'BondScale', p.Results.BondScale, ...
    'BondGraph', A, ...
    'SortMolecules', p.Results.SortMolecules);

ax = p.Results.Axes;
if isempty(ax)
    fig = figure;
    ax = axes(fig);
else
    fig = ancestor(ax, 'figure');
end

holdState = ishold(ax);
hold(ax, 'on');

siteHandles = gobjects(numel(molecules), 1);
bondHandles = cell(numel(molecules), 1);

for m = 1:numel(molecules)
    mol = molecules{m};
    X = mol.cart;

    labels = local_to_cell_column(mol.labels);
    [~, ~, C] = unique(labels, 'stable');

    siteHandles(m) = scatter3(ax, X(:,1), X(:,2), X(:,3), ...
        p.Results.MarkerSize, C, 'filled', ...
        'MarkerEdgeColor', [0 0 0]);

    if p.Results.ShowBonds
        bondHandles{m} = local_plot_molecule_bonds_angstrom(ax, X, labels, p.Results.BondScale);
    end

    com = mean(X, 1);
    text(ax, com(1), com(2), com(3), sprintf('  mol %d', m), ...
        'FontSize', 9, 'Color', [0 0 0]);
end

boxHandle = gobjects(0);
if p.Results.DrawBox
    boxHandle = viz.draw_cell_box(ax, S.lattice, ...
        'VectorsAs', 'rows', ...
        'Color', [0 0 0], ...
        'LineWidth', 1.0);
end

axis(ax, 'equal');
grid(ax, 'on');
xlabel(ax, 'x / Angstrom');
ylabel(ax, 'y / Angstrom');
zlabel(ax, 'z / Angstrom');
view(ax, 3);

ttl = char(string(p.Results.Title));
if isempty(ttl)
    ttl = 'Unwrapped unit-cell molecules';
end
title(ax, ttl);

if ~holdState
    hold(ax, 'off');
end

out = struct();
out.fig = fig;
out.ax = ax;
out.structure = S;
out.molecules = molecules;
out.bond_graph = A;
out.site_handles = siteHandles;
out.bond_handles = bondHandles;
out.box_handle = boxHandle;

end

% =========================================================================
% Helpers
% =========================================================================

function h = local_plot_molecule_bonds_angstrom(ax, X, labels, bondScale)

n = size(X, 1);
h = gobjects(0);
count = 0;

for i = 1:(n-1)
    for j = (i+1):n
        d = norm(X(j,:) - X(i,:));

        if io.bond_graph_tools.is_bonded(labels{i}, labels{j}, d, bondScale)
            count = count + 1;
            pts = [X(i,:); X(j,:)];
            h(count,1) = plot3(ax, pts(:,1), pts(:,2), pts(:,3), ...
                '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.0); %#ok<AGROW>
        end
    end
end

end

function c = local_to_cell_column(x)

if isstring(x)
    c = cellstr(x(:));
elseif iscellstr(x)
    c = x(:);
elseif iscell(x)
    c = x(:);
else
    error('viz:plot_unwrapped_unitcell_molecules:BadTextMetadata', ...
        'Text metadata must be string/cellstr/cell.');
end

end