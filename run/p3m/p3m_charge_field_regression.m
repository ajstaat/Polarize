%% run_p3m_charge_field_regression
% Compare P3M periodic charge external field against direct periodic Ewald.
%
% First milestone:
%
%   E_q^P3M ~= E_q^Ewald
%
% for fixed periodic charges.

clear; clc; close all;

fprintf('=== P3M charge-field regression ===\n');

%% ------------------------------------------------------------------------
% User controls

cfg = struct();

cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename   = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

cfg.supercellSize = [2 5 1];
cfg.bondScale     = 1.20;
cfg.pairCharges   = [+1 -1];

cfg.use_thole = false;  % keep first P3M regression pure point-charge Ewald
cfg.verbose = true;

cfg.ewald = struct();
cfg.ewald.alpha = 0.35;
cfg.ewald.rcut = 11.5;
cfg.ewald.kcut = 3.5;
cfg.ewald.boundary = 'tinfoil';

% Try anisotropic meshes roughly following cell shape.
meshList = { ...
    [24 40 24], ...
    [32 48 32], ...
    [40 64 40]};

orderList = [2 3 4];

%% ------------------------------------------------------------------------
% Import / build

fprintf('\n[setup] importing crystal template...\n');

crystal = io.import_contcar_as_crystal(cfg.filename, ...
    'BondScale', cfg.bondScale, ...
    'WrapFractional', true, ...
    'SortMolecules', false);

model = struct();
model.polarizable_classes = { ...
    'H_on_C_deg3', ...
    'H_on_C_deg4', ...
    'C_deg3', ...
    'C_deg4', ...
    'N', ...
    'O'};

model.alpha_by_class = struct( ...
    'H_on_C_deg3', 0.496, ...
    'H_on_C_deg4', 0.696, ...
    'C_deg3',      1.334, ...
    'C_deg4',      1.750, ...
    'N',           1.073, ...
    'O',           0.837);

model.thole_a = 0.39;

fprintf('\n[setup] building crystal system...\n');

buildOpts = struct();
buildOpts.supercell_size = cfg.supercellSize;
buildOpts.bondScale      = cfg.bondScale;
buildOpts.verbose        = cfg.verbose;
buildOpts.removeMolIDs   = [];
buildOpts.activeMolIDs   = [];

sys0 = builder.make_crystal_system(crystal, model, buildOpts);

fprintf('\n[setup] choosing charged-pair molecules...\n');

[refID, ~] = builder.choose_center_reference_molecule(sys0, ...
    'RequireComplete', true, ...
    'Verbose', cfg.verbose);

desc = builder.complete_molecule_descriptors_relative_to_reference(sys0, refID, ...
    'Verbose', cfg.verbose);

sameStack = builder.select_same_stack_neighbor(desc, 1, ...
    'Direction', 'either', ...
    'Verbose', cfg.verbose);

if isempty(sameStack.match_table) || height(sameStack.match_table) < 1
    error('Automatic same-stack shell-1 neighbor selection failed.');
end

nbrID = sameStack.match_table.molecule_id(1);

fprintf('\nChosen pair:\n');
fprintf('  ref ID = %d\n', refID);
fprintf('  nbr ID = %d\n', nbrID);

fprintf('\n[setup] applying charges and extracting polsys...\n');

sys = builder.apply_molecule_charges(sys0, [refID nbrID], ...
    'Mode', 'uniform', ...
    'TotalCharges', cfg.pairCharges, ...
    'SetActive', true, ...
    'DisablePolarizabilityOnCharged', true, ...
    'ZeroExistingCharges', true, ...
    'Verbose', cfg.verbose);

params_extract = struct();
params_extract.ewald = struct();
params_extract.ewald.mode = 'periodic';

polsys = builder.extract_polarization_system(sys, params_extract);
io.assert_atomic_units(polsys);

targetMask = logical(polsys.site_is_polarizable(:));
sourceMask = abs(polsys.site_charge(:)) > 0;

fprintf('\nSystem summary:\n');
fprintf('  n_sites       = %d\n', polsys.n_sites);
fprintf('  n_targets     = %d\n', nnz(targetMask));
fprintf('  n_sources     = %d\n', nnz(sourceMask));
fprintf('  total q source = %+0.16e\n', sum(polsys.site_charge(sourceMask)));

%% ------------------------------------------------------------------------
% Reference Ewald charge field

fprintf('\n[reference] building periodic Ewald charge field...\n');

fieldRef = struct();
fieldRef.mode = 'periodic';
fieldRef.exclude_self = true;
fieldRef.use_thole_damping = cfg.use_thole;
fieldRef.target_mask = targetMask;
fieldRef.source_mask = sourceMask;
fieldRef.verbose = true;
fieldRef.kspace_mode = 'full';
fieldRef.k_block_size = 2048;
fieldRef.kspace_memory_limit_gb = 8;
fieldRef.ewald = cfg.ewald;

tRef = tic;
[Eref, partsRef] = thole.induced_field_from_charges_periodic(polsys, fieldRef);
timeRef = toc(tRef);

fprintf('\nReference done:\n');
fprintf('  time        = %.6f s\n', timeRef);
fprintf('  ||Eref||_F  = %.16e\n', norm(Eref, 'fro'));
fprintf('  ||Ereal||_F = %.16e\n', norm(partsRef.real, 'fro'));
fprintf('  ||Erecip||_F= %.16e\n', norm(partsRef.recip, 'fro'));
fprintf('  ||Esurf||_F = %.16e\n', norm(partsRef.surf, 'fro'));

%% ------------------------------------------------------------------------
% P3M sweeps

rows = {};

fprintf('\n============================================================\n');
fprintf('P3M SWEEP\n');
fprintf('============================================================\n');

for im = 1:numel(meshList)
    meshSize = meshList{im};

    for io = 1:numel(orderList)
        order = orderList(io);

        fprintf('\n------------------------------------------------------------\n');
        fprintf('mesh = [%d %d %d], order = %d\n', meshSize, order);
        fprintf('------------------------------------------------------------\n');

        p3mOpts = struct();
        p3mOpts = struct();
        p3mOpts.ewald = cfg.ewald;
        p3mOpts.mesh_size = meshSize;
        p3mOpts.assignment_order = order;
        p3mOpts.target_mask = targetMask;
        p3mOpts.source_mask = sourceMask;
        p3mOpts.exclude_self = true;
        p3mOpts.use_thole_damping = cfg.use_thole;
        
        p3mOpts.realspace_backend = 'thole_periodic_real';
        p3mOpts.k_block_size = 2048;
        p3mOpts.kspace_memory_limit_gb = 8;
        
        p3mOpts.deconvolve_assignment = true;
        p3mOpts.deconvolution_floor = 1e-6;
        p3mOpts.verbose = true;

        tP3M = tic;
        [Ep3m, partsP3M] = p3m.compute_external_field_charges(polsys, p3mOpts);
        timeP3M = toc(tP3M);

        dE = Ep3m - Eref;

        targetSites = find(targetMask);
        errMag = vecnorm(dE(targetSites, :), 2, 2);
        refMag = vecnorm(Eref(targetSites, :), 2, 2);

        normRef = norm(Eref, 'fro');
        normP3M = norm(Ep3m, 'fro');
        normErr = norm(dE, 'fro');
        relErr = normErr / max(1e-300, normRef);

        rmsErr = sqrt(mean(errMag.^2));
        maxErr = max(errMag);
        rmsRef = sqrt(mean(refMag.^2));
        maxRef = max(refMag);

        fprintf('\nP3M result:\n');
        fprintf('  time          = %.6f s\n', timeP3M);
        fprintf('  ||Ep3m||_F    = %.16e\n', normP3M);
        fprintf('  ||Eref||_F    = %.16e\n', normRef);
        fprintf('  ||diff||_F    = %.16e\n', normErr);
        fprintf('  rel diff      = %.16e\n', relErr);
        fprintf('  target RMS err= %.16e\n', rmsErr);
        fprintf('  target max err= %.16e\n', maxErr);
        fprintf('  target RMS ref= %.16e\n', rmsRef);
        fprintf('  target max ref= %.16e\n', maxRef);

        rows(end+1,:) = { ...
            meshSize(1), meshSize(2), meshSize(3), order, ...
            normRef, normP3M, normErr, relErr, ...
            rmsRef, rmsErr, maxRef, maxErr, ...
            norm(partsP3M.real,'fro'), norm(partsP3M.recip,'fro'), norm(partsP3M.surf,'fro'), ...
            partsP3M.rho_total, partsP3M.solveInfo.nK, ...
            timeRef, timeP3M, partsP3M.time_real, partsP3M.time_mesh}; %#ok<SAGROW>
    end
end

summary = cell2table(rows, 'VariableNames', { ...
    'M1','M2','M3','order', ...
    'normRef','normP3M','normErr','relErr', ...
    'rmsRef','rmsErr','maxRef','maxErr', ...
    'normRealP3M','normRecipP3M','normSurfP3M', ...
    'rhoTotal','nKmesh', ...
    'timeRef','timeP3M','timeRealP3M','timeMeshP3M'});

fprintf('\n============================================================\n');
fprintf('P3M REGRESSION SUMMARY\n');
fprintf('============================================================\n');
disp(summary);

fprintf('\nDone.\n');