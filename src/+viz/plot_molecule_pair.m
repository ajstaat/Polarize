function out = plot_molecule_pair(sys, molIDs, varargin)
%PLOT_MOLECULE_PAIR Plot one molecule pair from a built system.
%
% out = viz.plot_molecule_pair(sys, [refMolID nbrMolID], Name, Value)
%
% Options
%   'Axes'          axes handle, default creates one
%   'DrawBox'       logical, default true
%   'ShowBonds'     logical, default true
%   'ShowCOM'       logical, default true
%   'RequireComplete' logical, default true
%   'Title'         plot title, default ''

p = inputParser;
addRequired(p, 'sys', @isstruct);
addRequired(p, 'molIDs', @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'Axes', [], @(x) isempty(x) || isgraphics(x, 'axes'));
addParameter(p, 'DrawBox', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ShowBonds', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ShowCOM', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'RequireComplete', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));
parse(p, sys, molIDs, varargin{:});

molIDs = reshape(molIDs, 1, 2);

if p.Results.RequireComplete
    local_assert_complete(sys, molIDs);
end

ax = p.Results.Axes;
if isempty(ax)
    fig = figure;
    ax = axes(fig);
else
    fig = ancestor(ax, 'figure');
end

holdState = ishold(ax);
hold(ax, 'on');

idx1 = builder.site_indices_for_molecule(sys, molIDs(1));
idx2 = builder.site_indices_for_molecule(sys, molIDs(2));

if isempty(idx1) || isempty(idx2)
    error('viz:plot_molecule_pair:MissingMolecule', ...
        'One or both molecule IDs were not found.');
end

out1 = local_scatter_molecule(ax, sys, idx1, 70);
out2 = local_scatter_molecule(ax, sys, idx2, 70);

bondHandle = gobjects(0);
if p.Results.ShowBonds
    bondHandle = [
        local_plot_molecule_bonds(ax, sys, idx1)
        local_plot_molecule_bonds(ax, sys, idx2)
    ];
end

boxHandle = gobjects(0);
if p.Results.DrawBox
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        boxHandle = viz.draw_cell_box(ax, sys.super_lattice, ...
            'VectorsAs', 'rows', ...
            'Color', [0 0 0], ...
            'LineWidth', 1.0);
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        boxHandle = viz.draw_cell_box(ax, sys.lattice, ...
            'VectorsAs', 'rows', ...
            'Color', [0 0 0], ...
            'LineWidth', 1.0);
    end
end

comHandle = gobjects(0);
pairLine = gobjects(0);

if p.Results.ShowCOM
    com1 = mean(sys.site_pos(idx1, :), 1);
    com2 = mean(sys.site_pos(idx2, :), 1);

    comHandle = scatter3(ax, [com1(1); com2(1)], [com1(2); com2(2)], [com1(3); com2(3)], ...
        100, 'kx', 'LineWidth', 2);

    pairLine = plot3(ax, [com1(1) com2(1)], [com1(2) com2(2)], [com1(3) com2(3)], ...
        'k--', 'LineWidth', 1.25);
end

axis(ax, 'equal');
grid(ax, 'on');
xlabel(ax, 'x / bohr');
ylabel(ax, 'y / bohr');
zlabel(ax, 'z / bohr');
view(ax, 3);

ttl = char(string(p.Results.Title));
if isempty(ttl)
    ttl = sprintf('Molecule pair %d - %d', molIDs(1), molIDs(2));
end
title(ax, ttl);

if ~holdState
    hold(ax, 'off');
end

out = struct();
out.fig = fig;
out.ax = ax;
out.mol_ids = molIDs;
out.site_indices = {idx1(:), idx2(:)};
out.site_handles = [out1; out2];
out.bond_handle = bondHandle;
out.box_handle = boxHandle;
out.com_handle = comHandle;
out.pair_line = pairLine;

end

% =========================================================================
% Helpers
% =========================================================================

function local_assert_complete(sys, molIDs)

if ~isfield(sys, 'molecule_table') || ...
        ~isfield(sys.molecule_table, 'molecule_id') || ...
        ~isfield(sys.molecule_table, 'is_complete_in_display')
    error('viz:plot_molecule_pair:MissingMoleculeTable', ...
        'sys.molecule_table with molecule_id and is_complete_in_display is required.');
end

T = sys.molecule_table;

for k = 1:numel(molIDs)
    row = find(T.molecule_id == molIDs(k), 1, 'first');

    if isempty(row)
        error('viz:plot_molecule_pair:UnknownMolecule', ...
            'Molecule ID %d was not found.', molIDs(k));
    end

    if ~T.is_complete_in_display(row)
        error('viz:plot_molecule_pair:IncompleteMolecule', ...
            'Molecule ID %d is not complete in the displayed supercell.', molIDs(k));
    end
end

end

function h = local_scatter_molecule(ax, sys, idx, markerSize)

X = sys.site_pos(idx, :);
types = local_to_cell_column(sys.site_type(idx));
[~, ~, C] = unique(types, 'stable');

h = scatter3(ax, X(:,1), X(:,2), X(:,3), markerSize, C, 'filled');
h.MarkerEdgeColor = [0 0 0];

end

function h = local_plot_molecule_bonds(ax, sys, idx)

X = sys.site_pos(idx, :);
types = local_to_cell_column(sys.site_type(idx));

toAng = local_length_to_angstrom_factor(sys);

n = numel(idx);
h = gobjects(0);
count = 0;

for i = 1:(n-1)
    for j = (i+1):n
        dAng = norm(X(j,:) - X(i,:)) * toAng;

        if io.bond_graph_tools.is_bonded(types{i}, types{j}, dAng, 1.20)
            count = count + 1;

            pts = [X(i,:); X(j,:)];
            h(count,1) = plot3(ax, pts(:,1), pts(:,2), pts(:,3), ...
                '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.0); %#ok<AGROW>
        end
    end
end

end

function factor = local_length_to_angstrom_factor(sys)

BOHR2ANG = 1 / 1.8897259886;

if isfield(sys, 'units') && isfield(sys.units, 'length')
    unit = lower(char(string(sys.units.length)));
else
    unit = 'bohr';
end

switch unit
    case {'angstrom', 'angstroms', 'ang', 'a'}
        factor = 1.0;

    case {'bohr', 'bohrs', 'a0', 'au_length'}
        factor = BOHR2ANG;

    otherwise
        error('viz:plot_molecule_pair:UnsupportedLengthUnit', ...
            'Unsupported length unit: %s', unit);
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
    error('viz:plot_molecule_pair:BadTextMetadata', ...
        'Text metadata must be string/cellstr/cell.');
end

end