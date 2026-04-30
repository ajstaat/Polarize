%% validate_aromatic_molecular_polarizabilities
%
% Validate molecular polarizabilities for benzene/naphthalene/anthracene
% and the actual 72-atom EP-PDI monomer from the VASP crystal structure.
%
% The acenes use idealized planar geometries.
%
% The EP-PDI monomer is extracted from the imported VASP crystal by building
% the molecular graph and selecting a 72-atom molecule.
%
% Units:
%   Coordinates: Angstrom internally converted to bohr.
%   Atomic polarizabilities: Angstrom^3 internally converted to bohr^3.
%   Field: atomic units.
%   Molecular polarizability output: Angstrom^3.

clear; clc; close all;

fprintf('\n============================================================\n');
fprintf('Validate aromatic / EP-PDI molecular polarizabilities\n');
fprintf('============================================================\n');

%% ------------------------------------------------------------------------
% Controls
% -------------------------------------------------------------------------

cfg = struct();

% Optional actual EP-PDI geometry from your crystal.
cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');
cfg.bondScale = 1.20;
cfg.ep_pdi_expected_atoms = 72;
cfg.ep_pdi_molecule_index = 1;   % if multiple 72-atom molecules are found

% Ideal acene geometry controls.
cfg.cc_bond_A = 1.397;
cfg.ch_bond_A = 1.083;

% AMOEBA/Thole controls.
cfg.thole_a = 0.39;

cfg.alpha_C_arom_A3 = 1.750;
cfg.alpha_H_arom_A3 = 0.696;

cfg.alpha_C_alkane_A3 = 1.334;
cfg.alpha_H_alkane_A3 = 0.496;

cfg.alpha_C_carbonyl_A3 = 1.334;
cfg.alpha_N_A3 = 1.073;
cfg.alpha_O_A3 = 0.837;

cfg.field_amplitude_au = 1e-4;
cfg.use_central_difference = true;

% Optional scalar test. Leave 1.0 for direct AMOEBA-like comparison.
cfg.alpha_scale = 1.0;

% Numerical solver settings.
cfg.solve_method = 'direct';
cfg.tol = 1e-12;
cfg.max_iter = 500;

fprintf('\nControls:\n');
fprintf('  VASP file             = %s\n', cfg.filename);
fprintf('  Thole a               = %.6f\n', cfg.thole_a);
fprintf('  alpha C aromatic      = %.6f Angstrom^3\n', cfg.alpha_C_arom_A3);
fprintf('  alpha H aromatic      = %.6f Angstrom^3\n', cfg.alpha_H_arom_A3);
fprintf('  alpha C alkane        = %.6f Angstrom^3\n', cfg.alpha_C_alkane_A3);
fprintf('  alpha H alkane        = %.6f Angstrom^3\n', cfg.alpha_H_alkane_A3);
fprintf('  alpha C carbonyl      = %.6f Angstrom^3\n', cfg.alpha_C_carbonyl_A3);
fprintf('  alpha N               = %.6f Angstrom^3\n', cfg.alpha_N_A3);
fprintf('  alpha O               = %.6f Angstrom^3\n', cfg.alpha_O_A3);
fprintf('  alpha scale           = %.6f\n', cfg.alpha_scale);
fprintf('  field amplitude       = %.6e a.u.\n', cfg.field_amplitude_au);
fprintf('  central difference    = %d\n', cfg.use_central_difference);

%% ------------------------------------------------------------------------
% Reference table for acenes
% -------------------------------------------------------------------------

ref = local_reference_table();

fprintf('\nReference AMOEBA acene values from Table 3 / Angstrom^3:\n');
disp(ref(:, {'molecule', 'AMOEBA_ax', 'AMOEBA_ay', 'AMOEBA_az'}));

%% ------------------------------------------------------------------------
% Build molecule list
% -------------------------------------------------------------------------

molecules = {};

molecules{end+1} = local_build_polyacene(1, cfg, "benzene");      %#ok<SAGROW>
molecules{end+1} = local_build_polyacene(2, cfg, "naphthalene");  %#ok<SAGROW>
molecules{end+1} = local_build_polyacene(3, cfg, "anthracene");   %#ok<SAGROW>

if isfile(cfg.filename)
    molecules{end+1} = local_build_ep_pdi_from_vasp(cfg);         %#ok<SAGROW>
else
    warning('VASP file not found. Skipping actual EP-PDI monomer:\n  %s', cfg.filename);
end

%% ------------------------------------------------------------------------
% Run molecules
% -------------------------------------------------------------------------

allSummary = table();

for m = 1:numel(molecules)
    mol = molecules{m};
    name = char(mol.name);

    fprintf('\n============================================================\n');
    fprintf('%s\n', upper(name));
    fprintf('============================================================\n');

    fprintf('  n atoms       = %d\n', mol.nAtoms);
    fprintf('  n carbon      = %d\n', nnz(strcmp(mol.element, 'C')));
    fprintf('  n hydrogen    = %d\n', nnz(strcmp(mol.element, 'H')));
    fprintf('  n nitrogen    = %d\n', nnz(strcmp(mol.element, 'N')));
    fprintf('  n oxygen      = %d\n', nnz(strcmp(mol.element, 'O')));

    if isfield(mol, 'type_label')
        fprintf('\n  Type counts:\n');
        local_print_type_counts(mol.type_label);
    end

    alphaTensor_A3 = local_compute_molecular_polarizability(mol, cfg);

    alphaSym_A3 = 0.5 .* (alphaTensor_A3 + alphaTensor_A3.');
    alphaAnti_A3 = 0.5 .* (alphaTensor_A3 - alphaTensor_A3.');

    [principalAxes, principalValsMat] = eig(alphaSym_A3);
    principalVals = diag(principalValsMat);

    [principalValsSorted, order] = sort(principalVals, 'descend');
    principalAxesSorted = principalAxes(:, order);

    computedDiag = [alphaSym_A3(1,1), alphaSym_A3(2,2), alphaSym_A3(3,3)];
    computedAvg = trace(alphaSym_A3) / 3;

    fprintf('\nComputed molecular polarizability tensor / Angstrom^3:\n');
    disp(alphaTensor_A3);

    fprintf('Symmetrized tensor / Angstrom^3:\n');
    disp(alphaSym_A3);

    fprintf('Antisymmetric norm ratio:\n');
    fprintf('  ||anti(alpha)||_F / ||sym(alpha)||_F = %.12e\n', ...
        norm(alphaAnti_A3, 'fro') / max(norm(alphaSym_A3, 'fro'), eps));

    fprintf('\nComputed diagonal in construction/crystal Cartesian axes [x y z] / Angstrom^3:\n');
    fprintf('  [%.6f %.6f %.6f]\n', computedDiag);

    fprintf('Computed alpha_avg / Angstrom^3:\n');
    fprintf('  %.6f\n', computedAvg);

    fprintf('Computed principal values sorted descending / Angstrom^3:\n');
    fprintf('  [%.6f %.6f %.6f]\n', principalValsSorted);

    fprintf('\nPrincipal axes columns:\n');
    disp(principalAxesSorted);

    summary = table();
    summary.molecule = string(name);
    summary.nAtoms = mol.nAtoms;
    summary.nC = nnz(strcmp(mol.element, 'C'));
    summary.nH = nnz(strcmp(mol.element, 'H'));
    summary.nN = nnz(strcmp(mol.element, 'N'));
    summary.nO = nnz(strcmp(mol.element, 'O'));

    summary.alpha_xx = alphaSym_A3(1,1);
    summary.alpha_yy = alphaSym_A3(2,2);
    summary.alpha_zz = alphaSym_A3(3,3);
    summary.alpha_xy = alphaSym_A3(1,2);
    summary.alpha_xz = alphaSym_A3(1,3);
    summary.alpha_yz = alphaSym_A3(2,3);
    summary.alpha_avg = computedAvg;

    summary.principal_1 = principalValsSorted(1);
    summary.principal_2 = principalValsSorted(2);
    summary.principal_3 = principalValsSorted(3);

    if any(strcmp(ref.molecule, string(name)))
        refRow = ref(strcmp(ref.molecule, string(name)), :);

        refAmoeba = [refRow.AMOEBA_ax, refRow.AMOEBA_ay, refRow.AMOEBA_az];

        summary.AMOEBA_ax = refAmoeba(1);
        summary.AMOEBA_ay = refAmoeba(2);
        summary.AMOEBA_az = refAmoeba(3);

        summary.err_x = computedDiag(1) - refAmoeba(1);
        summary.err_y = computedDiag(2) - refAmoeba(2);
        summary.err_z = computedDiag(3) - refAmoeba(3);

        summary.relerr_x_percent = 100 .* summary.err_x ./ refAmoeba(1);
        summary.relerr_y_percent = 100 .* summary.err_y ./ refAmoeba(2);
        summary.relerr_z_percent = 100 .* summary.err_z ./ refAmoeba(3);

        fprintf('\nReference AMOEBA Table 3 [x y z] / Angstrom^3:\n');
        fprintf('  [%.6f %.6f %.6f]\n', refAmoeba);

        fprintf('Difference computed diag - AMOEBA / Angstrom^3:\n');
        fprintf('  [%+.6f %+.6f %+.6f]\n', computedDiag - refAmoeba);

        fprintf('Relative difference computed diag vs AMOEBA:\n');
        fprintf('  [%+.3f%% %+.3f%% %+.3f%%]\n', ...
            100 .* (computedDiag - refAmoeba) ./ refAmoeba);
    else
        summary.AMOEBA_ax = NaN;
        summary.AMOEBA_ay = NaN;
        summary.AMOEBA_az = NaN;
        summary.err_x = NaN;
        summary.err_y = NaN;
        summary.err_z = NaN;
        summary.relerr_x_percent = NaN;
        summary.relerr_y_percent = NaN;
        summary.relerr_z_percent = NaN;
    end

    allSummary = [allSummary; summary]; %#ok<AGROW>
end

%% ------------------------------------------------------------------------
% Summary
% -------------------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('Overall summary\n');
fprintf('============================================================\n');

disp(allSummary);

fprintf('\nWorkflow completed successfully.\n');

%% =========================================================================
% Local functions
% =========================================================================

function mol = local_build_polyacene(nRings, cfg, name)
s = cfg.cc_bond_A;

anglesDeg = [30 90 150 210 270 330];
hexLocal = s .* [cosd(anglesDeg(:)), sind(anglesDeg(:)), zeros(6,1)];

centerSpacing = sqrt(3) * s;

C = zeros(0, 3);

for r = 1:nRings
    center = [(r-1) * centerSpacing, 0, 0];
    verts = hexLocal + center;
    C = local_append_unique_points(C, verts, 1e-8);
end

C = C - mean(C, 1);

nC = size(C, 1);
Cgraph = false(nC, nC);

for i = 1:nC
    for j = i+1:nC
        d = norm(C(i,:) - C(j,:));
        if abs(d - s) < 1e-5
            Cgraph(i,j) = true;
            Cgraph(j,i) = true;
        end
    end
end

degreeC = sum(Cgraph, 2);

H = zeros(0, 3);

for i = 1:nC
    if degreeC(i) ~= 2
        continue;
    end

    nbr = find(Cgraph(i,:));
    nbrMean = mean(C(nbr, :), 1);

    outward = C(i,:) - nbrMean;
    outward(3) = 0;
    outward = outward ./ norm(outward);

    H(end+1, :) = C(i,:) + cfg.ch_bond_A .* outward; %#ok<AGROW>
end

coordsA = [C; H];
element = [repmat({'C'}, nC, 1); repmat({'H'}, size(H,1), 1)];

alphaA3 = zeros(size(coordsA,1), 1);
typeLabel = strings(size(coordsA,1), 1);

alphaA3(strcmp(element, 'C')) = cfg.alpha_scale .* cfg.alpha_C_arom_A3;
alphaA3(strcmp(element, 'H')) = cfg.alpha_scale .* cfg.alpha_H_arom_A3;

typeLabel(strcmp(element, 'C')) = "C_aromatic";
typeLabel(strcmp(element, 'H')) = "H_aromatic";

mol = struct();
mol.name = string(name);
mol.nAtoms = size(coordsA, 1);
mol.coords_A = coordsA;
mol.element = element;
mol.alpha_A3 = alphaA3;
mol.type_label = typeLabel;
mol.nRings = nRings;
mol.Cgraph = Cgraph;
end

function mol = local_build_ep_pdi_from_vasp(cfg)
fprintf('\nExtracting actual EP-PDI monomer from VASP file...\n');

crystal = io.import_contcar_as_crystal(cfg.filename, ...
    'BondScale', cfg.bondScale, ...
    'SortMolecules', false);

model = local_dummy_model_for_molecule_build(cfg);

buildOpts = struct();
buildOpts.supercell_size = [1 1 1];
buildOpts.bondScale = cfg.bondScale;
buildOpts.verbose = false;

sys = builder.make_crystal_system(crystal, model, buildOpts);

molIDs = unique(sys.site_mol_id(:));
molIDs = molIDs(molIDs > 0);

counts = zeros(numel(molIDs), 1);

for k = 1:numel(molIDs)
    counts(k) = nnz(sys.site_mol_id == molIDs(k));
end

candidate = molIDs(counts == cfg.ep_pdi_expected_atoms);

if isempty(candidate)
    fprintf('Available molecule sizes in [1 1 1] system:\n');
    disp(table(molIDs(:), counts(:), 'VariableNames', {'molecule_id', 'nAtoms'}));

    error('Could not find a %d-atom molecule in the imported VASP system.', ...
        cfg.ep_pdi_expected_atoms);
end

if cfg.ep_pdi_molecule_index > numel(candidate)
    error('Requested EP-PDI molecule index %d, but only found %d candidates.', ...
        cfg.ep_pdi_molecule_index, numel(candidate));
end

molID = candidate(cfg.ep_pdi_molecule_index);
idx = find(sys.site_mol_id == molID);

coordsBohr = sys.site_pos(idx, :);
coordsA = coordsBohr ./ 1.88972612462577;

element = sys.site_type(idx);
element = element(:);

% Build a simple intramolecular bond graph from distances and covalent radii.
bondGraph = local_build_bond_graph_from_elements(coordsA, element, cfg.bondScale);

[alphaA3, typeLabel] = local_assign_ep_pdi_alpha(coordsA, element, bondGraph, cfg);

% Recenter for molecular polarizability. This does not change the response
% to a uniform field, but makes the tensor output easier to interpret.
coordsA = coordsA - mean(coordsA, 1);

mol = struct();
mol.name = "ep_pdi_vasp_monomer";
mol.nAtoms = numel(idx);
mol.coords_A = coordsA;
mol.element = element;
mol.alpha_A3 = alphaA3;
mol.type_label = typeLabel;
mol.bondGraph = bondGraph;
mol.source_molecule_id = molID;

fprintf('  selected molecule ID = %d\n', molID);
fprintf('  n atoms              = %d\n', mol.nAtoms);
end

function model = local_dummy_model_for_molecule_build(cfg)
% Minimal model only used to let builder.make_crystal_system identify
% molecules and carry site fields. These alpha values are overwritten for
% the extracted monomer.

model = struct();
model.thole_a = cfg.thole_a;
model.alpha_units = 'angstrom^3';

model.polarizable_classes = { ...
    'H_on_C_deg3', ...
    'H_on_C_deg4', ...
    'C_deg3', ...
    'C_deg4', ...
    'N', ...
    'O'};

model.alpha_by_class = struct( ...
    'H_on_C_deg3', cfg.alpha_H_arom_A3, ...
    'H_on_C_deg4', cfg.alpha_H_alkane_A3, ...
    'C_deg3', cfg.alpha_C_arom_A3, ...
    'C_deg4', cfg.alpha_C_alkane_A3, ...
    'N', cfg.alpha_N_A3, ...
    'O', cfg.alpha_O_A3);
end

function bondGraph = local_build_bond_graph_from_elements(coordsA, element, bondScale)
n = size(coordsA, 1);
bondGraph = false(n, n);

for i = 1:n
    for j = i+1:n
        ri = local_covalent_radius_A(element{i});
        rj = local_covalent_radius_A(element{j});

        cutoff = bondScale * (ri + rj);
        d = norm(coordsA(i,:) - coordsA(j,:));

        if d <= cutoff
            bondGraph(i,j) = true;
            bondGraph(j,i) = true;
        end
    end
end
end

function r = local_covalent_radius_A(el)
switch upper(strtrim(el))
    case 'H'
        r = 0.31;
    case 'C'
        r = 0.76;
    case 'N'
        r = 0.71;
    case 'O'
        r = 0.66;
    case 'S'
        r = 1.05;
    otherwise
        error('No covalent radius available for element "%s".', el);
end
end

function [alphaA3, typeLabel] = local_assign_ep_pdi_alpha(coordsA, element, bondGraph, cfg)
n = numel(element);

alphaA3 = zeros(n, 1);
typeLabel = strings(n, 1);

degree = sum(bondGraph, 2);

% First identify obvious carbonyl carbons: carbon bonded to at least one O.
isCarbonylC = false(n, 1);

for i = 1:n
    if ~strcmp(element{i}, 'C')
        continue;
    end

    nbr = find(bondGraph(i,:));
    if any(strcmp(element(nbr), 'O'))
        isCarbonylC(i) = true;
    end
end

% Carbon aromatic/sp2 heuristic:
%   - carbonyl C handled separately
%   - degree 3 carbon not carbonyl -> aromatic/sp2
%   - degree 4 carbon -> alkane/linker
%
% This is appropriate for the PDI core + phenyl rings + ethyl linkers.
isAromaticC = false(n, 1);
isAlkaneC = false(n, 1);

for i = 1:n
    if ~strcmp(element{i}, 'C')
        continue;
    end

    if isCarbonylC(i)
        continue;
    end

    if degree(i) <= 3
        isAromaticC(i) = true;
    else
        isAlkaneC(i) = true;
    end
end

% Hydrogens inherit aromatic/aliphatic type from their bonded heavy atom.
isAromaticH = false(n, 1);
isAlkaneH = false(n, 1);

for i = 1:n
    if ~strcmp(element{i}, 'H')
        continue;
    end

    nbr = find(bondGraph(i,:));

    if isempty(nbr)
        warning('Hydrogen atom %d has no bonded neighbor; assigning alkane H.', i);
        isAlkaneH(i) = true;
        continue;
    end

    parent = nbr(1);

    if isAromaticC(parent)
        isAromaticH(i) = true;
    else
        isAlkaneH(i) = true;
    end
end

for i = 1:n
    switch element{i}
        case 'C'
            if isCarbonylC(i)
                alphaA3(i) = cfg.alpha_scale .* cfg.alpha_C_carbonyl_A3;
                typeLabel(i) = "C_carbonyl";
            elseif isAromaticC(i)
                alphaA3(i) = cfg.alpha_scale .* cfg.alpha_C_arom_A3;
                typeLabel(i) = "C_aromatic";
            elseif isAlkaneC(i)
                alphaA3(i) = cfg.alpha_scale .* cfg.alpha_C_alkane_A3;
                typeLabel(i) = "C_alkane";
            else
                alphaA3(i) = cfg.alpha_scale .* cfg.alpha_C_alkane_A3;
                typeLabel(i) = "C_other_as_alkane";
            end

        case 'H'
            if isAromaticH(i)
                alphaA3(i) = cfg.alpha_scale .* cfg.alpha_H_arom_A3;
                typeLabel(i) = "H_aromatic";
            else
                alphaA3(i) = cfg.alpha_scale .* cfg.alpha_H_alkane_A3;
                typeLabel(i) = "H_alkane";
            end

        case 'N'
            alphaA3(i) = cfg.alpha_scale .* cfg.alpha_N_A3;
            typeLabel(i) = "N_imide";

        case 'O'
            alphaA3(i) = cfg.alpha_scale .* cfg.alpha_O_A3;
            typeLabel(i) = "O_carbonyl";

        otherwise
            error('Unsupported element "%s" in EP-PDI monomer.', element{i});
    end
end

if any(alphaA3 <= 0)
    bad = find(alphaA3 <= 0);
    error('Failed to assign positive polarizability to atoms: %s', mat2str(bad(:).'));
end

% Print useful typing diagnostics.
fprintf('\nEP-PDI polarizability typing diagnostics:\n');
local_print_type_counts(typeLabel);

fprintf('  total bare atomic alpha sum = %.6f Angstrom^3\n', sum(alphaA3));
fprintf('  degree counts:\n');
disp(tabulate(degree));
end

function local_print_type_counts(typeLabel)
u = unique(typeLabel);

for k = 1:numel(u)
    fprintf('    %-22s : %d\n', u(k), nnz(typeLabel == u(k)));
end
end

function P = local_append_unique_points(P, Q, tol)
for i = 1:size(Q,1)
    q = Q(i,:);

    if isempty(P)
        P = q;
        continue;
    end

    d = vecnorm(P - q, 2, 2);

    if all(d > tol)
        P(end+1,:) = q; %#ok<AGROW>
    end
end
end

function alphaMol_A3 = local_compute_molecular_polarizability(mol, cfg)
ANG3_TO_BOHR3 = 6.74833449394997;
BOHR_PER_ANG = 1.88972612462577;

coordsBohr = mol.coords_A .* BOHR_PER_ANG;
alphaBohr3 = mol.alpha_A3 .* ANG3_TO_BOHR3;

n = mol.nAtoms;
E0 = cfg.field_amplitude_au;

alphaMol_bohr3 = zeros(3, 3);

for j = 1:3
    Eplus = zeros(n, 3);
    Eplus(:, j) = E0;

    muPlus = local_solve_induced_dipoles(coordsBohr, alphaBohr3, Eplus, cfg);

    if cfg.use_central_difference
        Eminus = zeros(n, 3);
        Eminus(:, j) = -E0;

        muMinus = local_solve_induced_dipoles(coordsBohr, alphaBohr3, Eminus, cfg);

        Mresp = 0.5 .* (sum(muPlus, 1).' - sum(muMinus, 1).');
    else
        Mresp = sum(muPlus, 1).';
    end

    alphaMol_bohr3(:, j) = Mresp ./ E0;
end

alphaMol_A3 = alphaMol_bohr3 ./ ANG3_TO_BOHR3;
end

function mu = local_solve_induced_dipoles(coordsBohr, alphaBohr3, Eext, cfg)
n = size(coordsBohr, 1);
T = local_build_thole_dipole_tensor(coordsBohr, alphaBohr3, cfg.thole_a);

Avec = repelem(alphaBohr3(:), 3);
Evec = local_stack_xyz(Eext);

rhs = Avec .* Evec;
Aop = eye(3*n) - diag(Avec) * T;

switch lower(cfg.solve_method)
    case 'direct'
        muVec = Aop \ rhs;

    case 'gmres'
        [muVec, flag, relres] = gmres(Aop, rhs, [], cfg.tol, cfg.max_iter);

        if flag ~= 0
            warning('GMRES did not fully converge: flag=%d relres=%.3e', flag, relres);
        end

    otherwise
        error('Unsupported solve method "%s".', cfg.solve_method);
end

mu = local_unstack_xyz(muVec);
end

function T = local_build_thole_dipole_tensor(coordsBohr, alphaBohr3, tholeA)
n = size(coordsBohr, 1);
T = zeros(3*n, 3*n);

I3 = eye(3);

for i = 1:n
    for j = 1:n
        if i == j
            continue;
        end

        rvec = coordsBohr(i,:) - coordsBohr(j,:);
        r2 = dot(rvec, rvec);
        r = sqrt(r2);

        if r <= 0
            error('Coincident sites encountered.');
        end

        invR3 = 1 / (r^3);
        invR5 = 1 / (r^5);

        alphaIJ = alphaBohr3(i) * alphaBohr3(j);

        if alphaIJ > 0
            u = r / (alphaIJ^(1/6));
            damp = tholeA * u^3;

            expd = exp(-damp);
            f3 = 1 - expd;
            f5 = 1 - (1 + damp) * expd;
        else
            f3 = 1;
            f5 = 1;
        end

        block = 3 .* f5 .* (rvec(:) * rvec(:).') .* invR5 - ...
            f3 .* I3 .* invR3;

        rows = (3*(i-1)+1):(3*i);
        cols = (3*(j-1)+1):(3*j);

        T(rows, cols) = block;
    end
end
end

function v = local_stack_xyz(X)
n = size(X, 1);
v = reshape(X.', 3*n, 1);
end

function X = local_unstack_xyz(v)
v = v(:);
n = numel(v) / 3;

if n ~= round(n)
    error('Vector length must be divisible by 3.');
end

X = reshape(v, 3, n).';
end

function ref = local_reference_table()
ref = table();

ref.molecule = [
    "benzene"
    "naphthalene"
    "anthracene"
    ];

ref.Exp1_ax = [
    11.70
    20.20
    35.20
    ];

ref.Exp1_ay = [
    11.70
    18.80
    25.60
    ];

ref.Exp1_az = [
    5.72
    10.70
    15.20
    ];

ref.Exp2_ax = [
    12.26
    22.20
    44.70
    ];

ref.Exp2_ay = [
    12.26
    18.20
    25.80
    ];

ref.Exp2_az = [
    6.66
    7.30
    9.80
    ];

ref.AMOEBA_ax = [
    12.30
    21.78
    32.85
    ];

ref.AMOEBA_ay = [
    12.30
    18.51
    24.67
    ];

ref.AMOEBA_az = [
    6.64
    9.77
    12.63
    ];
end