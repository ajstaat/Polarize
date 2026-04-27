function h = draw_cell_box(ax, lattice, varargin)
%DRAW_CELL_BOX Draw a periodic cell box.
%
% h = viz.draw_cell_box(ax, lattice)
% h = viz.draw_cell_box(ax, lattice, Name, Value)
%
% Inputs
%   ax       target axes
%   lattice  3 x 3 direct lattice matrix
%
% Options
%   'Origin'     1 x 3 origin, default [0 0 0]
%   'VectorsAs'  'rows' or 'columns', default 'rows'
%   'Color'      line color, default [0 0 0]
%   'LineWidth'  default 1.0
%   'LineStyle'  default '-'
%
% Project convention:
%   lattice rows are direct vectors.

p = inputParser;
addRequired(p, 'ax');
addRequired(p, 'lattice', @(x) isnumeric(x) && isequal(size(x), [3 3]));
addParameter(p, 'Origin', [0 0 0], @(x) isnumeric(x) && numel(x) == 3);
addParameter(p, 'VectorsAs', 'rows', @(x) ischar(x) || isstring(x));
addParameter(p, 'Color', [0 0 0]);
addParameter(p, 'LineWidth', 1.0, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'LineStyle', '-', @(x) ischar(x) || isstring(x));
parse(p, ax, lattice, varargin{:});

origin = reshape(p.Results.Origin, 1, 3);
vectorsAs = lower(char(string(p.Results.VectorsAs)));

switch vectorsAs
    case 'rows'
        A = lattice(1,:);
        B = lattice(2,:);
        C = lattice(3,:);

    case {'columns', 'cols'}
        A = lattice(:,1).';
        B = lattice(:,2).';
        C = lattice(:,3).';

    otherwise
        error('viz:draw_cell_box:BadVectorsAs', ...
            'VectorsAs must be ''rows'' or ''columns''.');
end

corners = [
    origin
    origin + A
    origin + B
    origin + C
    origin + A + B
    origin + A + C
    origin + B + C
    origin + A + B + C
];

edges = [
    1 2
    1 3
    1 4
    2 5
    2 6
    3 5
    3 7
    4 6
    4 7
    5 8
    6 8
    7 8
];

holdState = ishold(ax);
hold(ax, 'on');

h = gobjects(size(edges,1), 1);

for k = 1:size(edges,1)
    pts = corners(edges(k,:), :);

    h(k) = plot3(ax, pts(:,1), pts(:,2), pts(:,3), ...
        'Color', p.Results.Color, ...
        'LineWidth', p.Results.LineWidth, ...
        'LineStyle', char(string(p.Results.LineStyle)));
end

if ~holdState
    hold(ax, 'off');
end

end