function test_viz_draw_cell_box()
%TEST_VIZ_DRAW_CELL_BOX Smoke test for viz.draw_cell_box.

fig = figure('Visible', 'off');
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>

ax = axes(fig);

H = [
    10 0 0
    1 11 0
    0 2 12
];

h = viz.draw_cell_box(ax, H, 'VectorsAs', 'rows');

assert(numel(h) == 12, ...
    'draw_cell_box should create 12 edge line handles.');

assert(all(isgraphics(h)), ...
    'draw_cell_box should return valid graphics handles.');

end