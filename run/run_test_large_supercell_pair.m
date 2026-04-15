clear; clc;

% ============================================================
% Test script:
%   1) build the generalized molecule-image descriptor table
%   2) select the nearest side-stack plane-parallel neighbor
%   3) plot the reference + that selected side-stack neighbor
%   4) move one step along that same stack by changing iy -> iy+1
%   5) plot the reference + the moved neighbor
%
% Assumption for this crystal:
%   - stack direction is along lattice axis b
%   - moving along the selected side stack means changing iy
%     while keeping base_mol_id, ix, and iz fixed
%
% Requires:
%   io.unwrap_all_contcar_molecules
%   builder.list_molecule_images_relative_to_reference
%   builder.select_molecule_image_neighbor
%   viz.plot_molecule_pair_from_contcar
% ============================================================

%% ------------------------------------------------------------
% File
% -------------------------------------------------------------
rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
filename = fullfile(rootFolder, 'a_0.0_CONTCAR.vasp');

%% ------------------------------------------------------------
% User settings
% -------------------------------------------------------------
bondScale = 1.20;
sortMolecules = false;

% Reference base molecule in the unit cell
refMolID = 1;
refShift = [0 0 0];

% Search / descriptor settings
searchImages = [4 4 4];
comCutoff = 20.0;

% For this crystal, stack direction is along b
axisMode = 'lattice';
axisSpec = 'b';

% Tolerances
sideStackPerpTol = 0.5;
planeAngleTolDeg = 10.0;

%% ------------------------------------------------------------
% Build generalized descriptor table
% -------------------------------------------------------------
fprintf('\n=== Building generalized molecule-image descriptor table ===\n');

Timg = builder.list_molecule_images_relative_to_reference(filename, ...
    'RefMolID', refMolID, ...
    'RefShift', refShift, ...
    'BondScale', bondScale, ...
    'SortMolecules', sortMolecules, ...
    'SearchImages', searchImages, ...
    'COMCutoff', comCutoff, ...
    'AxisMode', axisMode, ...
    'AxisSpec', axisSpec, ...
    'ComputePlaneNormals', true, ...
    'ComputeLongAxes', true);

disp(Timg(:, {'base_mol_id','ix','iy','iz','axis_parallel','axis_perp','plane_angle_deg','is_reference'}))

%% ------------------------------------------------------------
% Select nearest side-stack plane-parallel neighbor
% -------------------------------------------------------------
fprintf('\n=== Selecting nearest side-stack plane-parallel neighbor ===\n');

rowSideParallel = builder.select_molecule_image_neighbor(Timg, ...
    'Mode', 'nearest_side_stack', ...
    'BaseMolMode', 'different', ...
    'RequirePlaneParallel', true, ...
    'PlaneAngleTolDeg', planeAngleTolDeg, ...
    'AxisPerpTol', sideStackPerpTol);

disp('Selected side-stack plane-parallel anchor row:')
disp(rowSideParallel(:, {'base_mol_id','ix','iy','iz','axis_parallel','axis_perp','plane_angle_deg'}))

%% ------------------------------------------------------------
% Plot the selected side-stack anchor row
% -------------------------------------------------------------
fprintf('\n=== Plotting selected side-stack anchor row ===\n');

viz.plot_molecule_pair_from_contcar(filename, ...
    'RefMolID', refMolID, ...
    'RefShift', refShift, ...
    'NeighborRow', rowSideParallel, ...
    'BondScale', bondScale, ...
    'SortMolecules', sortMolecules, ...
    'CenterMode', 'midpoint', ...
    'ShowCOM', true, ...
    'ShowLabels', false, ...
    'ShowPairLine', true, ...
    'ShowTitle', true);

%% ------------------------------------------------------------
% Move one step along that same side stack:
% keep base_mol_id, ix, iz fixed; increment iy by +1
% -------------------------------------------------------------
fprintf('\n=== Moving one step along that same side stack (iy -> iy+1) ===\n');

rowMoved = Timg( ...
    Timg.base_mol_id == rowSideParallel.base_mol_id & ...
    Timg.ix == rowSideParallel.ix & ...
    Timg.iy == rowSideParallel.iy + 1 & ...
    Timg.iz == rowSideParallel.iz, :);

if isempty(rowMoved)
    error(['Could not find moved neighbor with same base_mol_id, ix, iz and iy+1. ' ...
           'Try increasing SearchImages or check whether this stack runs in the expected direction.']);
end

if height(rowMoved) > 1
    warning('Multiple rows matched the moved-neighbor query; using the first one.');
    rowMoved = rowMoved(1,:);
end

disp('Moved side-stack row:')
disp(rowMoved(:, {'base_mol_id','ix','iy','iz','axis_parallel','axis_perp','plane_angle_deg'}))

%% ------------------------------------------------------------
% Plot the moved side-stack row
% -------------------------------------------------------------
fprintf('\n=== Plotting moved side-stack row ===\n');

viz.plot_molecule_pair_from_contcar(filename, ...
    'RefMolID', refMolID, ...
    'RefShift', refShift, ...
    'NeighborRow', rowMoved, ...
    'BondScale', bondScale, ...
    'SortMolecules', sortMolecules, ...
    'CenterMode', 'midpoint', ...
    'ShowCOM', true, ...
    'ShowLabels', false, ...
    'ShowPairLine', true, ...
    'ShowTitle', true);

%% ------------------------------------------------------------
% Summary
% -------------------------------------------------------------
fprintf('\n=== Summary ===\n');

fprintf('Reference: b%d [%d %d %d]\n', ...
    refMolID, refShift(1), refShift(2), refShift(3));

fprintf('Selected side-stack anchor: b%d [%d %d %d], dpar=%.4f, dperp=%.4f, planeAngle=%.4f deg\n', ...
    rowSideParallel.base_mol_id, ...
    rowSideParallel.ix, rowSideParallel.iy, rowSideParallel.iz, ...
    rowSideParallel.axis_parallel, rowSideParallel.axis_perp, rowSideParallel.plane_angle_deg);

fprintf('Moved along stack:        b%d [%d %d %d], dpar=%.4f, dperp=%.4f, planeAngle=%.4f deg\n', ...
    rowMoved.base_mol_id, ...
    rowMoved.ix, rowMoved.iy, rowMoved.iz, ...
    rowMoved.axis_parallel, rowMoved.axis_perp, rowMoved.plane_angle_deg);