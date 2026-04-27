function test_viz_plot_unwrapped_unitcell_molecules()
%TEST_VIZ_PLOT_UNWRAPPED_UNITCELL_MOLECULES Smoke test for IO molecule plot.

[tmpFile, tmpDir] = local_write_tiny_poscar();
cleanup = onCleanup(@() local_cleanup(tmpDir)); %#ok<NASGU>

fig = figure('Visible', 'off');
cleanupFig = onCleanup(@() close(fig)); %#ok<NASGU>

ax = axes(fig);

out = viz.plot_unwrapped_unitcell_molecules(tmpFile, ...
    'Axes', ax, ...
    'ShowBonds', true, ...
    'DrawBox', true);

assert(isfield(out, 'molecules') && numel(out.molecules) == 1, ...
    'Expected one unwrapped molecule.');

assert(isfield(out, 'bond_graph') && nnz(triu(out.bond_graph,1)) == 1, ...
    'Expected one C-H bond.');

assert(isfield(out, 'ax') && isgraphics(out.ax), ...
    'Expected a valid axes handle.');

end

function [filename, tmpDir] = local_write_tiny_poscar()

tmpDir = tempname;
mkdir(tmpDir);

filename = fullfile(tmpDir, 'POSCAR_CH');

fid = fopen(filename, 'w');
assert(fid > 0, 'Failed to open temporary POSCAR.');

fprintf(fid, 'Tiny CH viz test\n');
fprintf(fid, '1.0\n');
fprintf(fid, '10.0 0.0 0.0\n');
fprintf(fid, '0.0 10.0 0.0\n');
fprintf(fid, '0.0 0.0 10.0\n');
fprintf(fid, 'C H\n');
fprintf(fid, '1 1\n');
fprintf(fid, 'Cartesian\n');
fprintf(fid, '5.00 5.00 5.00\n');
fprintf(fid, '6.09 5.00 5.00\n');

fclose(fid);

end

function local_cleanup(tmpDir)
if exist(tmpDir, 'dir')
    rmdir(tmpDir, 's');
end
end