%% run_vasp_tensor_clausius_mossotti_workflow
%
% Compute an anisotropic/tensor Clausius-Mossotti estimate of the crystal
% dielectric tensor from isolated molecular polarizability tensors.
%
% Workflow:
%   1. Import primitive VASP crystal.
%   2. Build [1 1 1] molecular system.
%   3. Extract each 72-atom EP-PDI molecule.
%   4. Assign AMOEBA-like atomic polarizabilities.
%   5. Compute isolated molecular alpha tensor in crystal Cartesian axes.
%   6. Sum alpha tensors over molecules in the primitive cell.
%   7. Apply tensor Clausius-Mossotti:
%
%        X = (4*pi / (3V)) * sum_m alpha_m
%
%        (eps - I) * inv(eps + 2I) = X
%
%      Therefore:
%
%        eps = (I + 2X) * inv(I - X)
%
% Units:
%   alpha tensors in bohr^3 for CM formula.
%   printed alpha also in Angstrom^3.
%   V in bohr^3.
%
% This is a local-field / Lorentz-Lorenz estimate, not a periodic SCF
% uniform-field calculation.

clear; clc; close all;

fprintf('\n============================================================\n');
fprintf('Tensor Clausius-Mossotti dielectric workflow\n');
fprintf('============================================================\n');

%% ------------------------------------------------------------------------
% Controls
% -------------------------------------------------------------------------

cfg = struct();

cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

cfg.bondScale = 1.20;
cfg.expected_atoms_per_molecule = 72;

% AMOEBA-like polarizabilities, Angstrom^3.
cfg.thole_a = 0.39;

cfg.alpha_C_arom_A3 = 1.750;
cfg.alpha_H_arom_A3 = 0.696;

cfg.alpha_C_alkane_A3 = 1.334;
cfg.alpha_H_alkane_A3 = 0.496;

cfg.alpha_C_carbonyl_A3 = 1.334;
cfg.alpha_N_A3 = 1.073;
cfg.alpha_O_A3 = 0.837;

cfg.alpha_scale = 1.0;

% Isolated molecule polarizability field amplitude.
cfg.field_amplitude_au = 1e-4;
cfg.use_central_difference = true;

% Solver for isolated molecular polarizabilities.
cfg.solve_method = 'direct';
cfg.tol = 1e-12;
cfg.max_iter = 500;

% Print per-molecule tensors.
cfg.print_molecule_details = true;

fprintf('\nInput structure:\n  %s\n', cfg.filename);
fprintf('\nControls:\n');
fprintf('  expected atoms/molecule = %d\n', cfg.expected_atoms_per_molecule);
fprintf('  Thole a                 = %.6f\n', cfg.thole_a);
fprintf('  alpha scale             = %.6f\n', cfg.alpha_scale);
fprintf('  field amplitude         = %.6e a.u.\n', cfg.field_amplitude_au);

if ~isfile(cfg.filename)
    error('Input file not found:\n  %s', cfg.filename);
end

%% ------------------------------------------------------------------------
% Constants
% -------------------------------------------------------------------------

ANG3_TO_BOHR3 = 6.74833449394997;
BOHR3_TO_ANG3 = 1.0 / ANG3_TO_BOHR3;

%% ------------------------------------------------------------------------
% 1. Import and build primitive molecular system
% -------------------------------------------------------------------------

fprintf('\n[1] Importing crystal template...\n');

crystal = io.import_contcar_as_crystal(cfg.filename, ...
    'BondScale', cfg.bondScale, ...
    'SortMolecules', false);

fprintf('  crystal.nSites    = %d\n', crystal.nSites);
fprintf('  crystal.nBaseMols = %d\n', crystal.nBaseMols);
fprintf('  crystal units     = %s\n', crystal.units.length);

fprintf('\n[2] Building [1 1 1] crystal system for molecule extraction...\n');

model = local_dummy_model_for_molecule_build(cfg);

buildOpts = struct();
buildOpts.supercell_size = [1 1 1];
buildOpts.bondScale = cfg.bondScale;
buildOpts.verbose = false;

sys = builder.make_crystal_system(crystal, model, buildOpts);
io.assert_atomic_units(sys);

lat = geom.get_lattice(sys);
V_bohr3 = lat.volume;
V_A3 = V_bohr3 * BOHR3_TO_ANG3;

fprintf('  sys.n_sites       = %d\n', sys.n_sites);
fprintf('  lattice volume    = %.12e bohr^3 = %.6f Angstrom^3\n', V_bohr3, V_A3);
fprintf('  H rows / bohr:\n');
disp(lat.H);

%% ------------------------------------------------------------------------
% 2. Extract molecules
% -------------------------------------------------------------------------

fprintf('\n[3] Extracting %d-atom molecules...\n', cfg.expected_atoms_per_molecule);

molIDs = unique(sys.site_mol_id(:));
molIDs = molIDs(molIDs > 0);

counts = zeros(numel(molIDs), 1);
for k = 1:numel(molIDs)
    counts(k) = nnz(sys.site_mol_id == molIDs(k));
end

candidateMolIDs = molIDs(counts == cfg.expected_atoms_per_molecule);

if isempty(candidateMolIDs)
    fprintf('\nAvailable molecule sizes:\n');
    disp(table(molIDs(:), counts(:), 'VariableNames', {'molecule_id', 'nAtoms'}));
    error('No %d-atom molecules found.', cfg.expected_atoms_per_molecule);
end

fprintf('  found %d candidate molecule(s): %s\n', ...
    numel(candidateMolIDs), mat2str(candidateMolIDs(:).'));

nMol = numel(candidateMolIDs);

%% ------------------------------------------------------------------------
% 3. Compute isolated molecular polarizability tensors
% -------------------------------------------------------------------------

fprintf('\n[4] Computing isolated molecular polarizability tensors...\n');

alphaMol_bohr3 = zeros(3, 3, nMol);
alphaMol_A3 = zeros(3, 3, nMol);

molSummary = table();

for im = 1:nMol
    molID = candidateMolIDs(im);
    idx = find(sys.site_mol_id == molID);

    mol = local_extract_molecule(sys, idx, cfg);

    [alpha_A3, alpha_bohr3] = local_compute_molecular_polarizability(mol, cfg);

    alphaMol_A3(:,:,im) = alpha_A3;
    alphaMol_bohr3(:,:,im) = alpha_bohr3;

    alphaSym_A3 = 0.5 .* (alpha_A3 + alpha_A3.');
    [eVec, eValMat] = eig(alphaSym_A3);
    vals = diag(eValMat);
    [valsSorted, order] = sort(vals, 'descend');
    eVec = eVec(:, order);

    if cfg.print_molecule_details
        fprintf('\nMolecule %d / ID %d:\n', im, molID);
        fprintf('  n atoms = %d\n', mol.nAtoms);
        fprintf('  type counts:\n');
        local_print_type_counts(mol.type_label);

        fprintf('  alpha tensor / Angstrom^3:\n');
        disp(alpha_A3);

        fprintf('  principal alpha / Angstrom^3 = [%.6f %.6f %.6f]\n', ...
            valsSorted(1), valsSorted(2), valsSorted(3));

        fprintf('  principal axes columns:\n');
        disp(eVec);
    end

    row = table();
    row.molecule_id = molID;
    row.nAtoms = mol.nAtoms;
    row.alpha_xx_A3 = alphaSym_A3(1,1);
    row.alpha_yy_A3 = alphaSym_A3(2,2);
    row.alpha_zz_A3 = alphaSym_A3(3,3);
    row.alpha_xy_A3 = alphaSym_A3(1,2);
    row.alpha_xz_A3 = alphaSym_A3(1,3);
    row.alpha_yz_A3 = alphaSym_A3(2,3);
    row.alpha_avg_A3 = trace(alphaSym_A3) / 3;
    row.principal_1_A3 = valsSorted(1);
    row.principal_2_A3 = valsSorted(2);
    row.principal_3_A3 = valsSorted(3);

    molSummary = [molSummary; row]; %#ok<AGROW>
end

fprintf('\nPer-molecule polarizability summary / Angstrom^3:\n');
disp(molSummary);

%% ------------------------------------------------------------------------
% 4. Tensor Clausius-Mossotti
% -------------------------------------------------------------------------

fprintf('\n[5] Computing tensor Clausius-Mossotti dielectric...\n');

alphaCell_bohr3 = sum(alphaMol_bohr3, 3);
alphaCell_A3 = alphaCell_bohr3 .* BOHR3_TO_ANG3;

alphaDensity = alphaCell_bohr3 ./ V_bohr3;

I3 = eye(3);
X = (4*pi/3) .* alphaDensity;

epsilonCM = (I3 + 2 .* X) / (I3 - X);
epsilonCMSym = 0.5 .* (epsilonCM + epsilonCM.');
epsilonCMAnti = 0.5 .* (epsilonCM - epsilonCM.');

[eVecCM, eValCM] = eig(epsilonCMSym);
epsPrincipal = diag(eValCM);
[epsPrincipalSorted, order] = sort(epsPrincipal, 'ascend');
eVecCM = eVecCM(:, order);

nPrincipal = sqrt(epsPrincipalSorted);

fprintf('\nCell polarizability tensor sum_m alpha_m / Angstrom^3:\n');
disp(alphaCell_A3);

fprintf('Cell polarizability tensor sum_m alpha_m / bohr^3:\n');
disp(alphaCell_bohr3);

fprintf('Polarizability density alpha_cell / V:\n');
disp(alphaDensity);

fprintf('X = (4*pi/3) alpha_cell/V:\n');
disp(X);

fprintf('\nTensor Clausius-Mossotti epsilon:\n');
disp(epsilonCM);

fprintf('Symmetrized tensor Clausius-Mossotti epsilon:\n');
disp(epsilonCMSym);

fprintf('Antisymmetric norm ratio:\n');
fprintf('  ||anti(eps)||_F / ||sym(eps)||_F = %.12e\n', ...
    norm(epsilonCMAnti, 'fro') / max(norm(epsilonCMSym, 'fro'), eps));

fprintf('\nPrincipal epsilon values, sorted ascending:\n');
fprintf('  [%.8f %.8f %.8f]\n', epsPrincipalSorted);

fprintf('Principal refractive indices sqrt(epsilon):\n');
fprintf('  [%.8f %.8f %.8f]\n', nPrincipal);

fprintf('Principal axes columns:\n');
disp(eVecCM);

%% ------------------------------------------------------------------------
% 5. Scalar Clausius-Mossotti sanity estimate
% -------------------------------------------------------------------------

fprintf('\n[6] Scalar isotropic Clausius-Mossotti sanity estimate...\n');

alphaAvgCell_bohr3 = trace(alphaCell_bohr3) / 3;
xScalar = (4*pi/3) * alphaAvgCell_bohr3 / V_bohr3;
epsScalar = (1 + 2*xScalar) / (1 - xScalar);
nScalar = sqrt(epsScalar);

fprintf('  alpha_cell_avg = %.6f Angstrom^3\n', alphaAvgCell_bohr3 * BOHR3_TO_ANG3);
fprintf('  alpha_mol_avg  = %.6f Angstrom^3 per molecule\n', ...
    alphaAvgCell_bohr3 * BOHR3_TO_ANG3 / nMol);
fprintf('  V per molecule = %.6f Angstrom^3\n', V_A3 / nMol);
fprintf('  x scalar       = %.12e\n', xScalar);
fprintf('  eps scalar CM  = %.8f\n', epsScalar);
fprintf('  n scalar CM    = %.8f\n', nScalar);

%% ------------------------------------------------------------------------
% 6. Compact summary
% -------------------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('Tensor Clausius-Mossotti summary\n');
fprintf('============================================================\n');

summary = table();
summary.nMolecules = nMol;
summary.volume_A3 = V_A3;
summary.volume_per_molecule_A3 = V_A3 / nMol;
summary.alpha_cell_avg_A3 = alphaAvgCell_bohr3 * BOHR3_TO_ANG3;
summary.alpha_mol_avg_A3 = alphaAvgCell_bohr3 * BOHR3_TO_ANG3 / nMol;
summary.eps_xx = epsilonCMSym(1,1);
summary.eps_yy = epsilonCMSym(2,2);
summary.eps_zz = epsilonCMSym(3,3);
summary.eps_xy = epsilonCMSym(1,2);
summary.eps_xz = epsilonCMSym(1,3);
summary.eps_yz = epsilonCMSym(2,3);
summary.eps_principal_1 = epsPrincipalSorted(1);
summary.eps_principal_2 = epsPrincipalSorted(2);
summary.eps_principal_3 = epsPrincipalSorted(3);
summary.n_principal_1 = nPrincipal(1);
summary.n_principal_2 = nPrincipal(2);
summary.n_principal_3 = nPrincipal(3);
summary.eps_scalar_CM = epsScalar;
summary.n_scalar_CM = nScalar;

disp(summary);

fprintf('\nWorkflow completed successfully.\n');

%% =========================================================================
% Local helpers
% =========================================================================

function model = local_dummy_model_for_molecule_build(cfg)
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

function mol = local_extract_molecule(sys, idx, cfg)
idx = idx(:);

coordsBohr = sys.site_pos(idx, :);
coordsA = coordsBohr ./ 1.88972612462577;

element = sys.site_type(idx);
element = element(:);

bondGraph = local_build_bond_graph_from_elements(coordsA, element, cfg.bondScale);

[alphaA3, typeLabel] = local_assign_ep_pdi_alpha(element, bondGraph, cfg);

% Recenter. Uniform-field polarizability is translation invariant.
coordsA = coordsA - mean(coordsA, 1);

mol = struct();
mol.nAtoms = numel(idx);
mol.coords_A = coordsA;
mol.element = element;
mol.alpha_A3 = alphaA3;
mol.type_label = typeLabel;
mol.bondGraph = bondGraph;
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

function [alphaA3, typeLabel] = local_assign_ep_pdi_alpha(element, bondGraph, cfg)
n = numel(element);

alphaA3 = zeros(n, 1);
typeLabel = strings(n, 1);

degree = sum(bondGraph, 2);

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
            error('Unsupported element "%s".', element{i});
    end
end

if any(alphaA3 <= 0)
    bad = find(alphaA3 <= 0);
    error('Failed to assign polarizability to atoms: %s', mat2str(bad(:).'));
end
end

function [alphaMol_A3, alphaMol_bohr3] = local_compute_molecular_polarizability(mol, cfg)
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

function local_print_type_counts(typeLabel)
u = unique(typeLabel);

for k = 1:numel(u)
    fprintf('    %-22s : %d\n', u(k), nnz(typeLabel == u(k)));
end
end