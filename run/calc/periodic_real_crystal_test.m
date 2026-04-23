%% test_periodic_external_field_real_crystal
% Focused test of calc.compute_external_field for a real charged-pair crystal.
%
% This script does NOT solve the induced-dipole SCF problem.
% It only builds nonperiodic and periodic external fields from fixed charges
% and checks:
%
%   1. mask/charge/alpha sanity
%   2. real vs reciprocal Ewald parts
%   3. alpha/kcut convergence behavior
%   4. full vs blocked k-space agreement
%   5. nonperiodic vs periodic qualitative scale
%
% Required recent edits:
%   - thole.induced_field_from_charges_periodic supports [E, parts]
%   - geom.build_target_source_field_cache_periodic uses
%     geom.shortest_lattice_translation(H)

clear; clc; close all;

fprintf('=== periodic external-field test: real charged-pair crystal ===\n');

%% ------------------------------------------------------------------------
% User controls

cfg = struct();

cfg.rootFolder = fullfile(getenv('HOME'), 'Desktop', 'Strain Spectra', 'structures');
cfg.filename   = fullfile(cfg.rootFolder, 'a_0.0_CONTCAR.vasp');

cfg.supercellSize = [2 5 1];
cfg.bondScale     = 1.20;

cfg.pairCharges = [+1 -1];

cfg.use_thole = true;
cfg.verbose = true;

% Real-space Ewald cutoff. Must obey rcut < Lmin/2 for the supercell.
cfg.rcut = 11.5;   % bohr

% A modest starting sweep. Increase kcut if field does not stabilize.
cfg.alphaList = [0.20 0.25 0.30 0.35 0.40 0.45 0.50];

% Option A: fixed kcut for simple first look.
cfg.kcut = 3.50;   % bohr^-1

% Option B: coupled alpha/kcut sweep, using kcut = kcutOverAlpha * alpha.
% Set cfg.useCoupledKcut = true to use this instead of fixed cfg.kcut.
cfg.useCoupledKcut = false;
cfg.kcutOverAlpha = 6.0;

cfg.boundary = 'tinfoil';

cfg.kspace_mode = 'full';
cfg.k_block_size = 2048;
cfg.kspace_memory_limit_gb = 8;

% Nonperiodic comparison. This intentionally uses uncut external field,
% matching the earlier nonperiodic example.
cfg.doNonperiodicComparison = true;

% Full-vs-blocked check can be expensive. It is useful on smaller kcut.
cfg.doFullVsBlockedCheck = false;
cfg.fullVsBlockedAlpha = 0.35;
cfg.fullVsBlockedKcut  = 1.00;

% Optional plots.
cfg.doPlots = true;
cfg.plotMaxTargets = 1000;

%% ------------------------------------------------------------------------
% Import crystal template

fprintf('\n[setup] importing crystal template...\n');

crystal = io.import_contcar_as_crystal(cfg.filename, ...
    'BondScale', cfg.bondScale, ...
    'WrapFractional', true, ...
    'SortMolecules', false);

fprintf('Imported crystal template:\n');
fprintf('  nSites     = %d\n', size(crystal.cart_coords, 1));
fprintf('  nBaseMols  = %d\n', numel(unique(crystal.base_mol_id)));

%% ------------------------------------------------------------------------
% Model

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

%% ------------------------------------------------------------------------
% Build working system

fprintf('\n[setup] building crystal system...\n');

buildOpts = struct();
buildOpts.supercell_size = cfg.supercellSize;
buildOpts.bondScale      = cfg.bondScale;
buildOpts.verbose        = cfg.verbose;
buildOpts.removeMolIDs   = [];
buildOpts.activeMolIDs   = [];

sys0 = builder.make_crystal_system(crystal, model, buildOpts);

fprintf('\nBuilt working system:\n');
fprintf('  n_unit_sites = %d\n', sys0.n_unit_sites);
fprintf('  n_cells      = %d\n', sys0.n_cells);
fprintf('  n_sites      = %d\n', sys0.n_sites);

%% ------------------------------------------------------------------------
% Automatically choose charged pair

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

%% ------------------------------------------------------------------------
% Apply charges and extract canonical polarization system

fprintf('\n[setup] applying charges and extracting polarization system...\n');

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
H = polsys.super_lattice;
pos = polsys.site_pos;

fprintf('units.length = %s\n', polsys.units.length);
fprintf('site_pos extent = [%g %g %g]\n', max(pos,[],1)-min(pos,[],1));
fprintf('H row lengths   = [%g %g %g]\n', norm(H(1,:)), norm(H(2,:)), norm(H(3,:)));
fprintf('H col lengths   = [%g %g %g]\n', norm(H(:,1)), norm(H(:,2)), norm(H(:,3)));

io.assert_atomic_units(polsys);

fprintf('\nCanonical polarization system:\n');
fprintf('  n_sites             = %d\n', polsys.n_sites);
fprintf('  n_polarizable_sites = %d\n', nnz(polsys.site_is_polarizable));

%% ------------------------------------------------------------------------
% Mask and alpha sanity

targetMask = logical(polsys.site_is_polarizable(:));
sourceMask = abs(polsys.site_charge(:)) > 0;

sourceCharge = polsys.site_charge(sourceMask);
sourceAlpha  = polsys.site_alpha(sourceMask);
targetAlpha  = polsys.site_alpha(targetMask);

fprintf('\n============================================================\n');
fprintf('MASK / CHARGE / ALPHA SANITY\n');
fprintf('============================================================\n');
fprintf('  n target sites (polarizable) = %d\n', nnz(targetMask));
fprintf('  n source sites (charged)     = %d\n', nnz(sourceMask));
fprintf('  target/source overlap        = %d\n', nnz(targetMask & sourceMask));
fprintf('  total source charge          = %+0.16e\n', sum(sourceCharge));
fprintf('  max |source charge|          = %.16e\n', max(abs(sourceCharge)));
fprintf('  min source alpha             = %.16e\n', min(sourceAlpha));
fprintf('  max source alpha             = %.16e\n', max(sourceAlpha));
fprintf('  min target alpha             = %.16e\n', min(targetAlpha));
fprintf('  max target alpha             = %.16e\n', max(targetAlpha));

if nnz(targetMask & sourceMask) ~= 0
    error('Charged source sites are still marked polarizable.');
end

if abs(sum(sourceCharge)) > 1e-10
    error('Selected charged sources are not neutral.');
end

if cfg.use_thole && any(sourceAlpha <= 0)
    error('use_thole = true, but at least one charged source site has nonpositive alpha.');
end

%% ------------------------------------------------------------------------
% Lattice / cutoff sanity

H = local_get_direct_lattice(polsys);
Lmin = geom.shortest_lattice_translation(H);
rcutMaxSafe = 0.5 * Lmin - 1e-12 * max(1, Lmin);

fprintf('\n============================================================\n');
fprintf('LATTICE / CUTOFF SANITY\n');
fprintf('============================================================\n');
fprintf('  supercell size    = [%d %d %d]\n', cfg.supercellSize);
fprintf('  volume            = %.16e bohr^3\n', abs(det(H)));
fprintf('  Lmin              = %.16e bohr\n', Lmin);
fprintf('  Lmin/2            = %.16e bohr\n', 0.5 * Lmin);
fprintf('  max safe rcut     = %.16e bohr\n', rcutMaxSafe);
fprintf('  requested rcut    = %.16e bohr\n', cfg.rcut);

if ~(cfg.rcut < rcutMaxSafe)
    error('cfg.rcut violates the single-image condition.');
end

%% ------------------------------------------------------------------------
% Nonperiodic comparison

E_nonper = [];
time_nonper = NaN;

if cfg.doNonperiodicComparison
    fprintf('\n============================================================\n');
    fprintf('NONPERIODIC EXTERNAL FIELD\n');
    fprintf('============================================================\n');

    params_nonper = struct();
    params_nonper.use_thole = cfg.use_thole;
    params_nonper.field = struct();
    params_nonper.field.mode = 'nonperiodic';
    params_nonper.field.exclude_self = true;
    params_nonper.field.use_thole_damping = cfg.use_thole;
    params_nonper.field.target_mask = targetMask;
    params_nonper.field.source_mask = sourceMask;

    t0 = tic;
    E_nonper = calc.compute_external_field(polsys, params_nonper);
    time_nonper = toc(t0);

    local_print_field_summary('nonperiodic', E_nonper, targetMask, time_nonper);
end

%% ------------------------------------------------------------------------
% Periodic alpha/kcut sweep

fprintf('\n============================================================\n');
fprintf('PERIODIC EWALD EXTERNAL FIELD SWEEP\n');
fprintf('============================================================\n');

nSweep = numel(cfg.alphaList);

sweep = table();
sweep.alpha = cfg.alphaList(:);
sweep.kcut = zeros(nSweep, 1);
sweep.nK = zeros(nSweep, 1);
sweep.nRealEntries = zeros(nSweep, 1);
sweep.normReal = zeros(nSweep, 1);
sweep.normRecip = zeros(nSweep, 1);
sweep.normTotal = zeros(nSweep, 1);
sweep.maxTargetField = zeros(nSweep, 1);
sweep.realOverTotal = zeros(nSweep, 1);
sweep.recipOverTotal = zeros(nSweep, 1);
sweep.relToLast = NaN(nSweep, 1);
sweep.timeReal = zeros(nSweep, 1);
sweep.timeRecip = zeros(nSweep, 1);
sweep.timeTotal = zeros(nSweep, 1);
sweep.storageMode = strings(nSweep, 1);

E_periodic_all = cell(nSweep, 1);
parts_all = cell(nSweep, 1);

for ia = 1:nSweep
    alpha = cfg.alphaList(ia);

    if cfg.useCoupledKcut
        kcut = cfg.kcutOverAlpha * alpha;
    else
        kcut = cfg.kcut;
    end

    params_per = local_make_periodic_field_params( ...
        cfg, targetMask, sourceMask, alpha, cfg.rcut, kcut, cfg.kspace_mode);

    fprintf('\n[periodic sweep %d/%d]\n', ia, nSweep);
    fprintf('  alpha = %.6f\n', alpha);
    fprintf('  rcut  = %.6f bohr\n', cfg.rcut);
    fprintf('  kcut  = %.6f bohr^-1\n', kcut);

    t0 = tic;
    E_per = calc.compute_external_field(polsys, params_per);
    time_wrapper = toc(t0);

    % Direct parts call, same fieldParams, so we can inspect real/recip.
    % This recomputes the field. Keeping it separate avoids changing the
    % calc.compute_external_field wrapper right now.
    [E_parts_total, parts] = thole.induced_field_from_charges_periodic( ...
        polsys, params_per.field);

    errParts = norm(E_per - E_parts_total, 'fro') / max(1, norm(E_parts_total, 'fro'));
    if errParts > 1e-12
        warning('Wrapper field and direct parts field differ: relerr = %.3e', errParts);
    end

    E_periodic_all{ia} = E_per;
    parts_all{ia} = parts;

    sweep.kcut(ia) = kcut;
    sweep.nK(ia) = parts.nK;
    sweep.nRealEntries(ia) = parts.nRealEntries;
    sweep.normReal(ia) = norm(parts.real, 'fro');
    sweep.normRecip(ia) = norm(parts.recip, 'fro');
    sweep.normTotal(ia) = norm(E_per, 'fro');
    sweep.maxTargetField(ia) = max(vecnorm(E_per(targetMask, :), 2, 2));
    sweep.realOverTotal(ia) = sweep.normReal(ia) / max(eps, sweep.normTotal(ia));
    sweep.recipOverTotal(ia) = sweep.normRecip(ia) / max(eps, sweep.normTotal(ia));
    sweep.timeReal(ia) = parts.time_real;
    sweep.timeRecip(ia) = parts.time_recip;
    sweep.timeTotal(ia) = time_wrapper;
    sweep.storageMode(ia) = string(parts.storage_mode);

    local_print_field_summary('periodic total', E_per, targetMask, time_wrapper);
    fprintf('  ||Ereal||_F       = %.16e\n', sweep.normReal(ia));
    fprintf('  ||Erecip||_F      = %.16e\n', sweep.normRecip(ia));
    fprintf('  real/total        = %.16e\n', sweep.realOverTotal(ia));
    fprintf('  recip/total       = %.16e\n', sweep.recipOverTotal(ia));
    fprintf('  nK                = %d\n', parts.nK);
    fprintf('  nRealEntries      = %d\n', parts.nRealEntries);
    fprintf('  storage mode      = %s\n', parts.storage_mode);
    fprintf('  time real         = %.6f s\n', parts.time_real);
    fprintf('  time recip        = %.6f s\n', parts.time_recip);
end

E_ref = E_periodic_all{end};
for ia = 1:nSweep
    sweep.relToLast(ia) = norm(E_periodic_all{ia} - E_ref, 'fro') / ...
        max(1e-300, norm(E_ref, 'fro'));
end

fprintf('\n============================================================\n');
fprintf('SWEEP SUMMARY\n');
fprintf('============================================================\n');
disp(sweep);

%% ------------------------------------------------------------------------
% Compare nonperiodic to final periodic field

if cfg.doNonperiodicComparison
    fprintf('\n============================================================\n');
    fprintf('NONPERIODIC VS FINAL PERIODIC\n');
    fprintf('============================================================\n');

    E_per_final = E_periodic_all{end};
    dE = E_per_final - E_nonper;

    fprintf('  ||E_nonper||_F      = %.16e\n', norm(E_nonper, 'fro'));
    fprintf('  ||E_periodic||_F    = %.16e\n', norm(E_per_final, 'fro'));
    fprintf('  ||difference||_F    = %.16e\n', norm(dE, 'fro'));
    fprintf('  rel difference      = %.16e\n', ...
        norm(dE, 'fro') / max(1e-300, norm(E_per_final, 'fro')));

    fprintf('  target RMS nonper   = %.16e\n', local_rms_vec(E_nonper(targetMask,:)));
    fprintf('  target RMS periodic = %.16e\n', local_rms_vec(E_per_final(targetMask,:)));
    fprintf('  target RMS diff     = %.16e\n', local_rms_vec(dE(targetMask,:)));
end

%% ------------------------------------------------------------------------
% Full-vs-blocked k-space agreement check

if cfg.doFullVsBlockedCheck
    fprintf('\n============================================================\n');
    fprintf('FULL VS BLOCKED K-SPACE CHECK\n');
    fprintf('============================================================\n');

    alpha = cfg.fullVsBlockedAlpha;
    kcut = cfg.fullVsBlockedKcut;

    params_full = local_make_periodic_field_params( ...
        cfg, targetMask, sourceMask, alpha, cfg.rcut, kcut, 'full');

    params_blocked = local_make_periodic_field_params( ...
        cfg, targetMask, sourceMask, alpha, cfg.rcut, kcut, 'blocked');

    fprintf('  alpha = %.6f\n', alpha);
    fprintf('  rcut  = %.6f bohr\n', cfg.rcut);
    fprintf('  kcut  = %.6f bohr^-1\n', kcut);

    [E_full, parts_full] = thole.induced_field_from_charges_periodic( ...
        polsys, params_full.field);

    [E_blocked, parts_blocked] = thole.induced_field_from_charges_periodic( ...
        polsys, params_blocked.field);

    relFullBlocked = norm(E_full - E_blocked, 'fro') / max(1e-300, norm(E_full, 'fro'));

    fprintf('  nK full       = %d\n', parts_full.nK);
    fprintf('  nK blocked    = %d\n', parts_blocked.nK);
    fprintf('  ||E_full||    = %.16e\n', norm(E_full, 'fro'));
    fprintf('  ||E_blocked|| = %.16e\n', norm(E_blocked, 'fro'));
    fprintf('  rel diff      = %.16e\n', relFullBlocked);

    if relFullBlocked > 1e-12
        warning('Full and blocked k-space paths differ by more than 1e-12.');
    end
end

%% ------------------------------------------------------------------------
% Optional plots

if cfg.doPlots
    fprintf('\n============================================================\n');
    fprintf('PLOTS\n');
    fprintf('============================================================\n');

    E_final = E_periodic_all{end};
    parts_final = parts_all{end};

    local_plot_field_magnitudes(E_final, parts_final, targetMask);

    if cfg.doNonperiodicComparison
        local_plot_periodic_vs_nonperiodic(E_nonper, E_final, targetMask);
    end

    local_plot_field_vectors(polsys, E_final, targetMask, sourceMask, cfg.plotMaxTargets);
end

fprintf('\nDone.\n');

%% ========================================================================
% Local helpers
%% ========================================================================

function params = local_make_periodic_field_params(cfg, targetMask, sourceMask, alpha, rcut, kcut, kspaceMode)
    params = struct();
    params.use_thole = cfg.use_thole;

    params.field = struct();
    params.field.mode = 'periodic';
    params.field.exclude_self = true;
    params.field.use_thole_damping = cfg.use_thole;
    params.field.target_mask = targetMask;
    params.field.source_mask = sourceMask;
    params.field.verbose = false;

    params.field.kspace_mode = kspaceMode;
    params.field.k_block_size = cfg.k_block_size;
    params.field.kspace_memory_limit_gb = cfg.kspace_memory_limit_gb;

    params.field.ewald = struct();
    params.field.ewald.alpha = alpha;
    params.field.ewald.rcut = rcut;
    params.field.ewald.kcut = kcut;
    params.field.ewald.boundary = cfg.boundary;
end

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('Missing direct lattice on system.');
    end
end

function local_print_field_summary(label, E, targetMask, elapsed)
    Et = E(targetMask, :);
    mag = vecnorm(Et, 2, 2);

    fprintf('  %s:\n', label);
    fprintf('    time             = %.6f s\n', elapsed);
    fprintf('    ||E||_F          = %.16e\n', norm(E, 'fro'));
    fprintf('    target RMS |E|   = %.16e\n', sqrt(mean(mag.^2)));
    fprintf('    target mean |E|  = %.16e\n', mean(mag));
    fprintf('    target max |E|   = %.16e\n', max(mag));
end

function y = local_rms_vec(X)
    mag = vecnorm(X, 2, 2);
    y = sqrt(mean(mag.^2));
end

function local_plot_field_magnitudes(E, parts, targetMask)
    Etot = E(targetMask, :);
    Ereal = parts.real(targetMask, :);
    Erecip = parts.recip(targetMask, :);

    mtot = vecnorm(Etot, 2, 2);
    mreal = vecnorm(Ereal, 2, 2);
    mrecip = vecnorm(Erecip, 2, 2);

    figure('Name', 'Periodic external field magnitudes');
    tiledlayout(3,1);

    nexttile;
    histogram(mreal, 80);
    xlabel('|Ereal|');
    ylabel('count');
    title('Real-space field magnitude');

    nexttile;
    histogram(mrecip, 80);
    xlabel('|Erecip|');
    ylabel('count');
    title('Reciprocal-space field magnitude');

    nexttile;
    histogram(mtot, 80);
    xlabel('|Etotal|');
    ylabel('count');
    title('Total periodic field magnitude');
end

function local_plot_periodic_vs_nonperiodic(E_nonper, E_per, targetMask)
    En = E_nonper(targetMask, :);
    Ep = E_per(targetMask, :);

    mn = vecnorm(En, 2, 2);
    mp = vecnorm(Ep, 2, 2);
    md = vecnorm(Ep - En, 2, 2);

    figure('Name', 'Nonperiodic vs periodic external field');
    tiledlayout(2,1);

    nexttile;
    plot(mn, mp, '.');
    xlabel('|E nonperiodic|');
    ylabel('|E periodic|');
    title('Field magnitude comparison');
    grid on;
    axis equal;

    nexttile;
    histogram(md, 80);
    xlabel('|E periodic - E nonperiodic|');
    ylabel('count');
    title('Difference magnitude on target sites');
end

function local_plot_field_vectors(sys, E, targetMask, sourceMask, maxTargets)
    pos = sys.site_pos;

    targetIdx = find(targetMask);
    sourceIdx = find(sourceMask);

    if numel(targetIdx) > maxTargets
        % Deterministic thinning.
        keep = round(linspace(1, numel(targetIdx), maxTargets));
        targetIdxPlot = targetIdx(keep);
    else
        targetIdxPlot = targetIdx;
    end

    figure('Name', 'Periodic external field vectors');
    hold on;

    scatter3(pos(targetIdxPlot,1), pos(targetIdxPlot,2), pos(targetIdxPlot,3), ...
        8, 'filled');

    scatter3(pos(sourceIdx,1), pos(sourceIdx,2), pos(sourceIdx,3), ...
        60, 'filled');

    Eplot = E(targetIdxPlot, :);
    mag = vecnorm(Eplot, 2, 2);
    nonzero = mag > 0;

    if any(nonzero)
        scale = 0.2 * local_box_scale(sys) / max(mag(nonzero));
        quiver3( ...
            pos(targetIdxPlot(nonzero),1), ...
            pos(targetIdxPlot(nonzero),2), ...
            pos(targetIdxPlot(nonzero),3), ...
            scale * Eplot(nonzero,1), ...
            scale * Eplot(nonzero,2), ...
            scale * Eplot(nonzero,3), ...
            0);
    end

    if exist('viz.draw_cell_box', 'file') == 2 || exist('+viz/draw_cell_box.m', 'file') == 2
        H = local_get_direct_lattice(sys);
        H = H * 1.8897259886;
        frac = geom.cart_to_frac(sys.site_pos, H);
        frac_wrapped = frac - floor(frac);
        pos = geom.frac_to_cart(frac_wrapped, H);
        boxOrigin = [0 0 0];
        viz.draw_cell_box(gca, H, boxOrigin, 'rows');
    end

    xlabel('x / bohr');
    ylabel('y / bohr');
    zlabel('z / bohr');
    title('Periodic external field vectors');
    axis equal;
    grid on;
    view(3);
    hold off;
end

function s = local_box_scale(sys)
    H = local_get_direct_lattice(sys);
    s = max([norm(H(1,:)), norm(H(2,:)), norm(H(3,:))]);
end