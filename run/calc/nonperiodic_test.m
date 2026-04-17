%% run_test_nonperiodic_uniform_charged_pair_calc
clear; clc; close all;

rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
filename   = fullfile(rootFolder, 'a_0.0_CONTCAR.vasp');

%% ------------------------------------------------------------------------
% Import crystal template
crystal = io.import_contcar_as_crystal(filename, ...
    'BondScale', 1.20, ...
    'WrapFractional', true, ...
    'SortMolecules', false);

fprintf('Imported crystal template:\n');
fprintf('  nSites     = %d\n', size(crystal.cart_coords, 1));
fprintf('  nBaseMols  = %d\n', numel(unique(crystal.base_mol_id)));

%% ------------------------------------------------------------------------
% Model with graph-based polarizability classes
% Updated values:
%   C_deg3      -> 1.334
%   H_on_C_deg3 -> 0.496
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

% Refactored run_polarization_calc interface for uniform molecular charges
model.total_charges = [+1; -1];

%% ------------------------------------------------------------------------
% Build-time / workflow options
opts = struct();
opts.supercell_size = [2 5 1];
opts.bondScale = 1.20;
opts.verbose = true;
opts.removeMolIDs = [];

% Same-stack shell = 1 pair from your existing test
refID = 32;
nbrID = 38;
opts.activeMolIDs = [refID nbrID];

% Charged molecules act as fixed field sources only
opts.depolarizeActiveMolecules = true;

fprintf('\nChosen pair:\n');
fprintf('  ref ID = %d\n', refID);
fprintf('  nbr ID = %d\n', nbrID);

%% ------------------------------------------------------------------------
% Calculation parameters
params = struct();

params.output = struct();
params.output.verbose = true;
params.output.computeEnergy = true;
params.output.computeEnergyByMolecule = false;
params.output.computeDipoleDipoleDecomp = false;

params.field = struct();
params.field.include_external_charges = true;
params.field.exclude_self = true;
params.field.softening = 0.0;

params.scf = struct();
params.scf.solver = 'direct';
params.scf.tol = 1e-10;
params.scf.maxIter = 500;
params.scf.mixing = 0.5;
params.scf.omega = 1.0;
params.scf.verbose = true;
params.scf.printEvery = 25;
params.scf.residualEvery = 25;
params.scf.stopMetric = 'max_dmu';
params.scf.use_thole = true;
params.scf.softening = 0.0;

params.ewald = struct();
params.ewald.mode = 'nonperiodic';

%% ------------------------------------------------------------------------
% Run full nonperiodic polarization calculation
result = calc.run_polarization_calc(crystal, model, opts, params);

sys    = result.sys;
polsys = result.polsys;
Eext   = result.Eext;
Tpol   = result.Tpol;
problem = result.problem;
mu     = result.mu;
energy = result.energy;

%% ------------------------------------------------------------------------
% Diagnostics on charged pair
refIdx = builder.site_indices_for_molecule(sys, refID);
nbrIdx = builder.site_indices_for_molecule(sys, nbrID);

fprintf('\nPost-run diagnostics:\n');
fprintf('  ref site count                  = %d\n', numel(refIdx));
fprintf('  nbr site count                  = %d\n', numel(nbrIdx));
fprintf('  total system charge             = %+0.10f\n', sum(sys.site_charge));
fprintf('  ref total charge                = %+0.10f\n', sum(sys.site_charge(refIdx)));
fprintf('  nbr total charge                = %+0.10f\n', sum(sys.site_charge(nbrIdx)));
fprintf('  ref polarizable sites remaining = %d\n', nnz(sys.site_is_polarizable(refIdx)));
fprintf('  nbr polarizable sites remaining = %d\n', nnz(sys.site_is_polarizable(nbrIdx)));
fprintf('  active sites total              = %d\n', nnz(sys.site_is_active));

%% ------------------------------------------------------------------------
% Assertions: charged pair setup
if abs(sum(sys.site_charge(refIdx)) - 1.0) > 1e-12
    error('Reference molecule charge assignment failed.');
end

if abs(sum(sys.site_charge(nbrIdx)) + 1.0) > 1e-12
    error('Neighbor molecule charge assignment failed.');
end

if any(sys.site_is_polarizable(refIdx))
    error('Reference molecule still has polarizable sites enabled.');
end

if any(sys.site_is_polarizable(nbrIdx))
    error('Neighbor molecule still has polarizable sites enabled.');
end

if ~all(sys.site_is_active(refIdx)) || ~all(sys.site_is_active(nbrIdx))
    error('Active mask was not set correctly on charged molecules.');
end

%% ------------------------------------------------------------------------
% Assertions: canonical polarization system
io.assert_atomic_units(polsys);

if ~strcmpi(polsys.units.length, 'bohr')
    error('polsys length units are not bohr.');
end
if ~strcmpi(polsys.units.alpha, 'atomic_unit')
    error('polsys alpha units are not atomic_unit.');
end
if ~strcmpi(polsys.units.charge, 'elementary_charge')
    error('polsys charge units are not elementary_charge.');
end

if size(polsys.site_pos, 1) ~= polsys.n_sites
    error('polsys.site_pos row count mismatch.');
end

if numel(polsys.site_alpha) ~= polsys.n_sites
    error('polsys.site_alpha length mismatch.');
end

if numel(polsys.site_charge) ~= polsys.n_sites
    error('polsys.site_charge length mismatch.');
end

%% ------------------------------------------------------------------------
% Assertions: external field / active-space setup
if isempty(Eext)
    error('Expected nonempty external field.');
end

if ~isequal(size(Eext), [polsys.n_sites, 3])
    error('Eext must be N x 3.');
end

if any(~isfinite(Eext(:)))
    error('Eext contains non-finite values.');
end

if problem.nSites ~= polsys.n_sites
    error('problem.nSites does not match polsys.n_sites.');
end

if problem.nPolSites <= 0
    error('Expected at least one polarizable site.');
end

if numel(problem.activeSites) ~= problem.nPolSites
    error('problem.activeSites length mismatch.');
end

if numel(problem.alpha_pol_vec) ~= 3 * problem.nPolSites
    error('problem.alpha_pol_vec length mismatch.');
end

if numel(problem.Eext_pol_vec) ~= 3 * problem.nPolSites
    error('problem.Eext_pol_vec length mismatch.');
end

%% ------------------------------------------------------------------------
% Assertions: operator / induced dipoles
if isempty(Tpol)
    error('Expected nonempty Tpol.');
end

if ~isequal(size(Tpol), [3*problem.nPolSites, 3*problem.nPolSites])
    error('Tpol size mismatch.');
end

if any(~isfinite(Tpol(:)))
    error('Tpol contains non-finite values.');
end

if ~isequal(size(mu), [polsys.n_sites, 3])
    error('mu must be N x 3.');
end

if any(~isfinite(mu(:)))
    error('mu contains non-finite values.');
end

if ~isfield(result, 'direct') || ~isstruct(result.direct)
    error('Missing direct-solver diagnostics.');
end

if ~isfield(result.direct, 'relres') || ~isfinite(result.direct.relres)
    error('Direct-solver relres is missing or non-finite.');
end

if result.direct.relres > 1e-8
    error('Direct-solver relres too large: %.3e', result.direct.relres);
end

relres_check = thole.compute_active_space_relres(problem, Tpol, mu);
if ~isfinite(relres_check) || relres_check > 1e-8
    error('Independent active-space residual too large: %.3e', relres_check);
end

%% ------------------------------------------------------------------------
% Assertions: energy
if ~isfield(energy, 'external_charge_dipole') || ~isfinite(energy.external_charge_dipole)
    error('energy.external_charge_dipole must be finite.');
end

if ~isfield(energy, 'dipole_dipole') || ~isfinite(energy.dipole_dipole)
    error('energy.dipole_dipole must be finite.');
end

if ~isfield(energy, 'total') || ~isfinite(energy.total)
    error('energy.total must be finite.');
end

energy_check = calc.compute_total_energy_active_space( ...
    polsys, problem, mu, Eext, Tpol);

if abs(energy.total - energy_check.total) > 1e-12
    error('Energy mismatch between stored and recomputed total.');
end

energy_total_eV = result.energy.total * 27.211386245988;

%% ------------------------------------------------------------------------
% Optional plot
fig = figure(1);
clf(fig);
fig.Color = 'w';
ax = axes('Parent', fig);

viz.plot_supercell_selection(sys, ...
    'ReferenceMolID', refID, ...
    'NeighborMolID', nbrID, ...
    'RequireCompleteRef', true, ...
    'RequireCompleteNbr', true, ...
    'IncompleteAction', 'error', ...
    'EnvMarkerSize', 52, ...
    'RefMarkerSize', 78, ...
    'NbrMarkerSize', 78, ...
    'ShowCOM', true, ...
    'ShowLabels', true, ...
    'DrawBox', true, ...
    'Axes', ax, ...
    'Title', sprintf('Nonperiodic charged pair: ref %d (+1), nbr %d (-1)', refID, nbrID));

%% ------------------------------------------------------------------------
% Report
fprintf('\nCheckpoint passed.\n');
fprintf('  polsys n_sites             = %d\n', polsys.n_sites);
fprintf('  active polarizable sites   = %d\n', problem.nPolSites);
fprintf('  total assigned charge      = %+0.10f\n', sum(polsys.site_charge));
fprintf('  direct relres              = %.3e\n', result.direct.relres);
fprintf('  total polarization energy  = %+0.10f\n', energy.total);
fprintf('  total polarization energy  = %+0.10f eV\n', energy_total_eV);
fprintf('run_nonperiodic_uniform_charged_pair_calc completed successfully.\n');