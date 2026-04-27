function out = plot_system(sys, varargin)
%PLOT_SYSTEM Plot a molecule-aware system.
%
% out = viz.plot_system(sys)
% out = viz.plot_system(sys, Name, Value)
%
% Options
%   'Axes'          axes handle, default creates one
%   'MoleculeIDs'   molecule IDs to plot, default all
%   'ColorBy'       'site_type'|'molecule'|'complete'|'charge', default 'site_type'
%   'MarkerSize'    default 28
%   'Alpha'         marker face alpha, default 0.35
%   'DrawBox'       logical, default true
%   'ShowCOM'       logical, default true
%   'ShowBonds'     logical, default false
%   'Title'         plot title, default ''
%
% Output fields include:
%   .fig
%   .ax
%   .site_handle
%   .box_handle
%   .com_handle
%   .bond_handle

p = inputParser;
addRequired(p, 'sys', @isstruct);
addParameter(p, 'Axes', [], @(x) isempty(x) || isgraphics(x, 'axes'));
addParameter(p, 'MoleculeIDs', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'ColorBy', 'site_type', @(x) ischar(x) || isstring(x));
addParameter(p, 'MarkerSize', 28, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'Alpha', 0.35, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'DrawBox', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ShowCOM', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ShowBonds', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));
parse(p, sys, varargin{:});

validate_sys(sys);

ax = p.Results.Axes;
if isempty(ax)
    fig = figure;
    ax = axes(fig);
else
    fig = ancestor(ax, 'figure');
end

molIDs = p.Results.MoleculeIDs;
if isempty(molIDs)
    siteMask = true(size(sys.site_pos,1), 1);
else
    siteMask = ismember(sys.site_mol_id(:), molIDs(:));
end

X = sys.site_pos(siteMask, :);
siteIdx = find(siteMask);

colorBy = lower(char(string(p.Results.ColorBy)));
C = local_color_values(sys, siteIdx, colorBy);

holdState = ishold(ax);
hold(ax, 'on');

siteHandle = scatter3(ax, X(:,1), X(:,2), X(:,3), ...
    p.Results.MarkerSize, C, 'filled');

try
    siteHandle.MarkerFaceAlpha = p.Results.Alpha;
    siteHandle.MarkerEdgeAlpha = min(1, p.Results.Alpha + 0.25);
catch
end

boxHandle = gobjects(0);
if p.Results.DrawBox && isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
    boxHandle = viz.draw_cell_box(ax, sys.super_lattice, ...
        'VectorsAs', 'rows', ...
        'Color', [0 0 0], ...
        'LineWidth', 1.0);
elseif p.Results.DrawBox && isfield(sys, 'lattice') && ~isempty(sys.lattice)
    boxHandle = viz.draw_cell_box(ax, sys.lattice, ...
        'VectorsAs', 'rows', ...
        'Color', [0 0 0], ...
        'LineWidth', 1.0);
end

comHandle = gobjects(0);
if p.Results.ShowCOM && isfield(sys, 'molecule_table') && ~isempty(sys.molecule_table)
    T = sys.molecule_table;

    if isempty(molIDs)
        rows = 1:numel(T.molecule_id);
    else
        rows = find(ismember(T.molecule_id(:), molIDs(:)));
    end

    if ~isempty(rows) && isfield(T, 'com')
        COM = T.com(rows, :);
        comHandle = scatter3(ax, COM(:,1), COM(:,2), COM(:,3), ...
            60, 'kx', 'LineWidth', 1.5);
    end
end

bondHandle = gobjects(0);
if p.Results.ShowBonds
    bondHandle = local_plot_bonds(ax, sys, siteIdx);
end

axis(ax, 'equal');
grid(ax, 'on');
xlabel(ax, 'x / bohr');
ylabel(ax, 'y / bohr');
zlabel(ax, 'z / bohr');

ttl = char(string(p.Results.Title));
if ~isempty(ttl)
    title(ax, ttl);
end

view(ax, 3);

if ~holdState
    hold(ax, 'off');
end

out = struct();
out.fig = fig;
out.ax = ax;
out.site_handle = siteHandle;
out.box_handle = boxHandle;
out.com_handle = comHandle;
out.bond_handle = bondHandle;
out.site_indices = siteIdx;

end

% =========================================================================
% Helpers
% =========================================================================

function validate_sys(sys)

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos) || size(sys.site_pos,2) ~= 3
    error('viz:plot_system:BadSitePos', ...
        'sys.site_pos must be N x 3.');
end

n = size(sys.site_pos, 1);

if ~isfield(sys, 'site_mol_id') || numel(sys.site_mol_id) ~= n
    error('viz:plot_system:MissingSiteMolID', ...
        'sys.site_mol_id must have one entry per site.');
end

if ~isfield(sys, 'site_type') || numel(sys.site_type) ~= n
    error('viz:plot_system:MissingSiteType', ...
        'sys.site_type must have one entry per site.');
end

end

function C = local_color_values(sys, siteIdx, colorBy)

switch colorBy
    case 'molecule'
        C = double(sys.site_mol_id(siteIdx));

    case 'complete'
        C = zeros(numel(siteIdx), 1);

        if isfield(sys, 'molecule_table') && ...
                isfield(sys.molecule_table, 'molecule_id') && ...
                isfield(sys.molecule_table, 'is_complete_in_display')
            T = sys.molecule_table;

            for k = 1:numel(siteIdx)
                molID = sys.site_mol_id(siteIdx(k));
                row = find(T.molecule_id == molID, 1, 'first');

                if ~isempty(row)
                    C(k) = double(T.is_complete_in_display(row));
                end
            end
        end

    case 'charge'
        if isfield(sys, 'site_charge') && numel(sys.site_charge) >= max(siteIdx)
            C = double(sys.site_charge(siteIdx));
        else
            C = zeros(numel(siteIdx), 1);
        end

    case {'site_type', 'element'}
        types = local_to_cell_column(sys.site_type(siteIdx));
        [~, ~, C] = unique(types, 'stable');

    otherwise
        error('viz:plot_system:BadColorBy', ...
            'Unsupported ColorBy value: %s', colorBy);
end

end

function h = local_plot_bonds(ax, sys, siteIdx)

X = sys.site_pos(siteIdx, :);
types = local_to_cell_column(sys.site_type(siteIdx));

toAng = local_length_to_angstrom_factor(sys);

n = numel(siteIdx);
h = gobjects(0);
count = 0;

for i = 1:(n-1)
    for j = (i+1):n
        if sys.site_mol_id(siteIdx(i)) ~= sys.site_mol_id(siteIdx(j))
            continue;
        end

        dAng = norm(X(j,:) - X(i,:)) * toAng;

        if io.bond_graph_tools.is_bonded(types{i}, types{j}, dAng, 1.20)
            count = count + 1;

            pts = [X(i,:); X(j,:)];
            h(count,1) = plot3(ax, pts(:,1), pts(:,2), pts(:,3), ...
                '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.75); %#ok<AGROW>
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
        error('viz:plot_system:UnsupportedLengthUnit', ...
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
    error('viz:plot_system:BadTextMetadata', ...
        'Text metadata must be string/cellstr/cell.');
end

end