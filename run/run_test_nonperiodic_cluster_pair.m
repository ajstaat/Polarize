clear; clc;

%% ------------------------------------------------------------
% File
% -------------------------------------------------------------
rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
filename = fullfile(rootFolder, 'a_0.0_CONTCAR.vasp');

%% ------------------------------------------------------------
% General settings
% -------------------------------------------------------------
bondScale = 1.20;
sortMolecules = false;

refMolID = 1;
refShift = [0 0 0];

searchImages = [4 4 4];
comCutoff = 20.0;

axisMode = 'lattice';
axisSpec = 'b';

sideStackPerpTol = 0.5;
planeAngleTolDeg = 10.0;

clusterCutoff = 18.0;

%% ------------------------------------------------------------
% Build generalized descriptor table
% -------------------------------------------------------------
fprintf('\n=== Building molecule-image descriptor table ===\n');

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

%% ------------------------------------------------------------
% Select a side-stack plane-parallel anchor
% -------------------------------------------------------------
rowSideParallel = builder.select_molecule_image_neighbor(Timg, ...
    'Mode', 'nearest_side_stack', ...
    'BaseMolMode', 'different', ...
    'RequirePlaneParallel', true, ...
    'PlaneAngleTolDeg', planeAngleTolDeg, ...
    'AxisPerpTol', sideStackPerpTol);

disp('Selected side-stack plane-parallel anchor row:')
disp(rowSideParallel(:, {'base_mol_id','ix','iy','iz','axis_parallel','axis_perp','plane_angle_deg'}))

%% ------------------------------------------------------------
% Move one step along that same side stack (iy -> iy+1)
% -------------------------------------------------------------
rowMoved = Timg( ...
    Timg.base_mol_id == rowSideParallel.base_mol_id & ...
    Timg.ix == rowSideParallel.ix & ...
    Timg.iy == rowSideParallel.iy + 1 & ...
    Timg.iz == rowSideParallel.iz, :);

if isempty(rowMoved)
    error('Could not find moved neighbor by stepping iy -> iy+1.');
end
if height(rowMoved) > 1
    rowMoved = rowMoved(1,:);
end

disp('Moved side-stack row:')
disp(rowMoved(:, {'base_mol_id','ix','iy','iz','axis_parallel','axis_perp','plane_angle_deg'}))

%% ------------------------------------------------------------
% Build a nonperiodic cluster around the selected pair
% -------------------------------------------------------------
fprintf('\n=== Building nonperiodic cluster around the selected pair ===\n');

cluster = builder.build_nonperiodic_cluster_from_pair(filename, refMolID, rowMoved, ...
    'RefShift', refShift, ...
    'BondScale', bondScale, ...
    'SortMolecules', sortMolecules, ...
    'SearchImages', searchImages, ...
    'ClusterCutoff', clusterCutoff, ...
    'CenterMode', 'pair_midpoint', ...
    'BoxPadding', 5.0);

disp(cluster.focus_pair)
disp(cluster.image_table(:, {'unique_mol_id','base_mol_id','ix','iy','iz','cx','cy','cz'}))

%% ------------------------------------------------------------
% Build model
% -------------------------------------------------------------
model = struct();
model.polarizable_types = {'H','C','N','O'};
model.alpha_by_type = struct( ...
    'H', 0.696, ...
    'C', 1.750, ...
    'N', 1.073, ...
    'O', 0.837);
model.thole_a = 0.39;

% Whole-molecule charge distributions over the two selected cluster molecules
refUID = cluster.focus_pair.ref_unique_mol_id;
nbrUID = cluster.focus_pair.neighbor_unique_mol_id;

idxRef = (cluster.crystal.mol_id == refUID);
idxNbr = (cluster.crystal.mol_id == nbrUID);

labelsRef = cluster.crystal.site_label(idxRef);
labelsNbr = cluster.crystal.site_label(idxNbr);

nRef = numel(labelsRef);
nNbr = numel(labelsNbr);

model.charge_patterns = {
    struct('mol_id', refUID, ...
           'site_label', {labelsRef(:).'}, ...
           'delta_q',    (-1/nRef) * ones(1, nRef)), ...
    struct('mol_id', nbrUID, ...
           'site_label', {labelsNbr(:).'}, ...
           'delta_q',    (+1/nNbr) * ones(1, nNbr))
};

%% ------------------------------------------------------------
% Nonperiodic calc options
% -------------------------------------------------------------
opts = struct();
opts.supercell_size = [1 1 1];   % irrelevant for fragment crystal
opts.verbose = true;
opts.removeMolIDs = [];
opts.activeMolIDs = [refUID, nbrUID];
opts.depolarizeActiveMolecules = true;

params = util.default_params();
params.scf.solver = 'gs';
params.scf.verbose = true;
params.scf.printEvery = 100;
params.ewald.mode = 'nonperiodic';
params.field.softening = 0.5;
params.scf.softening = 0.5;

%% ------------------------------------------------------------
% Run nonperiodic polarization calculation on the cluster
% -------------------------------------------------------------
fprintf('\n=== Running nonperiodic polarization calculation on cluster ===\n');

result = calc.run_polarization_calc(cluster.crystal, model, opts, params);

disp(result.energy)