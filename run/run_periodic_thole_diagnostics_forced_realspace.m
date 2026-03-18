clear; clc;

% ============================================================
% FORCE A GEOMETRY WHERE REAL-SPACE CONTRIBUTIONS EXIST
% ============================================================
%
% Strategy:
%   - Build a much smaller unit cell than before so some periodic images
%     and intermolecular separations fall inside the automatically chosen rcut
%   - Keep the same project workflow
%   - Use automatic Ewald parameter selection from choose_ewald_params
%   - Compare bare periodic vs real-space-Thole-corrected periodic operators
%
% Goal:
%   Confirm that the real-space Thole correction is alive when the
%   real-space Ewald channel is actually nonzero.

% -----------------------------
% Define a compact crystal
% -----------------------------
crystal = struct();

% Smaller triclinic cell so rcut from automatic selection will actually
% include some short real-space interactions.
crystal.cellpar = [6.5, 7.0, 6.2, 88.0, 101.0, 84.0];
crystal.lattice = [];

% Put two molecules with one polarizable-heavy atom each and one auxiliary
% site each, arranged so some intermolecular distances are short.
%
% Molecule 1: C1, H1
% Molecule 2: N1, H2
%
% Fractional positions chosen to create several ~3-5 Å separations after
% converting to Cartesian coordinates.
crystal.frac_coords = [
    0.10 0.15 0.20   % site 1: C1  (mol 1)
    0.16 0.18 0.24   % site 2: H1  (mol 1)
    0.42 0.20 0.28   % site 3: N1  (mol 2)
    0.48 0.24 0.31   % site 4: H2  (mol 2)
];

crystal.cart_coords = [];
crystal.mol_id = [1; 1; 2; 2];
crystal.site_label = {'C1'; 'H1'; 'N1'; 'H2'};
crystal.site_type  = {'C';  'H';  'N';  'H' };

% -----------------------------
% Define model
% -----------------------------
model = struct();
model.polarizable_types = {'C', 'N'};
model.alpha_by_type = struct();
model.alpha_by_type.C = 10.0;
model.alpha_by_type.N = 8.0;
model.thole_a = 0.39;

% Put charge pattern on active molecules
model.charge_pattern = struct();
model.charge_pattern.site_label = {'C1', 'H1'};
model.charge_pattern.delta_q = [0.10, -0.10];

% -----------------------------
% Build neutral system first for molecule lookup
% -----------------------------
opts0 = struct();
opts0.supercell_size = [2 2 1];
opts0.verbose = true;
opts0.removeMolIDs = [];
opts0.activeMolIDs = [];

params0 = util.default_params();
res0 = calc.run_polarization_calc(crystal, model, opts0, params0);
sys0 = res0.sys;

disp('Available molecules:')
disp(builder.list_molecules(sys0))

% Choose two copies of base molecule 1 as active
molA = builder.find_unique_molecule_id(sys0, 1, [0 0 0]);
molB = builder.find_unique_molecule_id(sys0, 1, [1 0 0]);

opts = struct();
opts.supercell_size = [2 2 1];
opts.verbose = true;
opts.removeMolIDs = [];
opts.activeMolIDs = [molA, molB];

% ============================================================
% BARE PERIODIC OPERATOR
% ============================================================
params_bare = util.default_params();
params_bare.ewald.mode = 'periodic_triclinic';
params_bare.ewald.auto = true;
params_bare.ewald.tol = 1e-8;
params_bare.ewald.boundary = 'tinfoil';
params_bare.ewald.rcut_fraction = 0.9;
params_bare.ewald.use_thole_real_space = false;
params_bare.scf.solver = 'direct';

res_bare = calc.run_polarization_calc(crystal, model, opts, params_bare);

% ============================================================
% THOLE-CORRECTED PERIODIC OPERATOR
% ============================================================
params_thole = util.default_params();
params_thole.ewald.mode = 'periodic_triclinic';
params_thole.ewald.auto = true;
params_thole.ewald.tol = 1e-8;
params_thole.ewald.boundary = 'tinfoil';
params_thole.ewald.rcut_fraction = 0.9;
params_thole.ewald.use_thole_real_space = true;
params_thole.ewald.thole_a = model.thole_a;
params_thole.scf.solver = 'direct';

res_thole = calc.run_polarization_calc(crystal, model, opts, params_thole);

% ============================================================
% REPORT
% ============================================================
fprintf('\n============================================================\n');
fprintf('EWALD PARAMETERS USED\n');
fprintf('============================================================\n');
fprintf('alpha = %.16e\n', res_bare.Tmeta.alpha);
fprintf('rcut  = %.16e\n', res_bare.Tmeta.rcut);
fprintf('kcut  = %.16e\n', res_bare.Tmeta.kcut);
fprintf('num_kvec = %d\n', res_bare.Tmeta.num_kvec);
fprintf('real image bounds = [%d %d %d]\n', res_bare.Tmeta.real_image_bounds);

fprintf('\n============================================================\n');
fprintf('TOTAL ENERGY / DIPOLES\n');
fprintf('============================================================\n');
fprintf('Bare total energy         = %+ .15e\n', res_bare.energy.total);
fprintf('Thole total energy        = %+ .15e\n', res_thole.energy.total);
fprintf('|delta U_total|           = %.3e\n', abs(res_thole.energy.total - res_bare.energy.total));

dmu = res_thole.mu - res_bare.mu;
site_dmu = sqrt(sum(dmu.^2, 2));

fprintf('max |delta mu_i|          = %.16e\n', max(site_dmu));
fprintf('||delta mu||_F            = %.16e\n', norm(dmu, 'fro'));

fprintf('\n============================================================\n');
fprintf('REAL / RECIP / TOTAL MATRIX CHANGES\n');
fprintf('============================================================\n');

Treal_bare  = res_bare.Tparts.real;
Treal_thole = res_thole.Tparts.real;
dTreal = Treal_thole - Treal_bare;

fprintf('||Treal_bare||_F          = %.16e\n', norm(Treal_bare, 'fro'));
fprintf('||Treal_thole||_F         = %.16e\n', norm(Treal_thole, 'fro'));
fprintf('||dTreal||_F              = %.16e\n', norm(dTreal, 'fro'));
if norm(Treal_bare, 'fro') > 0
    fprintf('relative real-space change = %.16e\n', norm(dTreal,'fro') / norm(Treal_bare,'fro'));
else
    fprintf('relative real-space change = NaN (bare norm is zero)\n');
end

fprintf('||Trecip_bare||_F         = %.16e\n', norm(res_bare.Tparts.recip, 'fro'));
fprintf('||Trecip_thole||_F        = %.16e\n', norm(res_thole.Tparts.recip, 'fro'));
fprintf('||dTrecip||_F             = %.16e\n', norm(res_thole.Tparts.recip - res_bare.Tparts.recip, 'fro'));

fprintf('||T_bare||_F              = %.16e\n', norm(res_bare.T, 'fro'));
fprintf('||T_thole||_F             = %.16e\n', norm(res_thole.T, 'fro'));
fprintf('||dT_total||_F            = %.16e\n', norm(res_thole.T - res_bare.T, 'fro'));
if norm(res_bare.T, 'fro') > 0
    fprintf('relative total change      = %.16e\n', ...
        norm(res_thole.T - res_bare.T, 'fro') / norm(res_bare.T, 'fro'));
else
    fprintf('relative total change      = NaN (bare norm is zero)\n');
end

fprintf('\n============================================================\n');
fprintf('DIPOLE-DIPOLE ENERGY DECOMPOSITION\n');
fprintf('============================================================\n');
disp('Bare dipole-dipole decomposition:')
disp(res_bare.dipole_dipole_decomposition)

disp('Thole dipole-dipole decomposition:')
disp(res_thole.dipole_dipole_decomposition)

fprintf('delta U_dd(active,active) = %+ .15e\n', ...
    res_thole.dipole_dipole_decomposition.active_active - ...
    res_bare.dipole_dipole_decomposition.active_active);

fprintf('delta U_dd(active,env)    = %+ .15e\n', ...
    res_thole.dipole_dipole_decomposition.active_environment - ...
    res_bare.dipole_dipole_decomposition.active_environment);

fprintf('delta U_dd(env,env)       = %+ .15e\n', ...
    res_thole.dipole_dipole_decomposition.environment_environment - ...
    res_bare.dipole_dipole_decomposition.environment_environment);

fprintf('\n============================================================\n');
fprintf('CLOSEST POLARIZABLE SITE PAIRS\n');
fprintf('============================================================\n');

sys = res_bare.sys;
pos = sys.site_pos;
alpha_site = sys.site_alpha(:);
polMask = logical(sys.site_is_polarizable(:));
site_labels = sys.site_label;
site_mol_id = sys.site_mol_id(:);

polIdx = find(polMask);
nPol = numel(polIdx);

pairs = [];
for a = 1:nPol
    i = polIdx(a);
    for b = a+1:nPol
        j = polIdx(b);

        rij = pos(i,:) - pos(j,:);
        dist = norm(rij);

        pairs = [pairs; i, j, dist]; %#ok<AGROW>
    end
end

if isempty(pairs)
    fprintf('No polarizable pairs found.\n');
else
    pairs = sortrows(pairs, 3);
    nShow = min(12, size(pairs,1));

    fprintf('Showing %d closest polarizable pairs:\n', nShow);
    fprintf('   i    j      dist(A)        mol_i   mol_j      alpha_i      alpha_j    label_i   label_j\n');

    for p = 1:nShow
        i = pairs(p,1);
        j = pairs(p,2);
        dist = pairs(p,3);

        fprintf('%4d %4d   %12.6f   %6d %6d   %10.6f   %10.6f   %-6s   %-6s\n', ...
            i, j, dist, site_mol_id(i), site_mol_id(j), alpha_site(i), alpha_site(j), ...
            site_labels{i}, site_labels{j});
    end
end

fprintf('\n============================================================\n');
fprintf('THOLE f3/f5 FACTORS FOR CLOSEST POLARIZABLE PAIRS\n');
fprintf('============================================================\n');

if ~isempty(pairs)
    nShow = min(12, size(pairs,1));
    fprintf('   i    j      dist(A)            s               f3               f5               l3               l5\n');

    for p = 1:nShow
        i = pairs(p,1);
        j = pairs(p,2);
        dist = pairs(p,3);

        tf = thole.thole_f3f5_factors(dist, alpha_site(i), alpha_site(j), model.thole_a);

        fprintf('%4d %4d   %12.6f   %14.6e   %14.6e   %14.6e   %14.6e   %14.6e\n', ...
            i, j, dist, tf.s, tf.f3, tf.f5, tf.l3, tf.l5);
    end
end

fprintf('\n============================================================\n');
fprintf('PAIRWISE BLOCK CHANGE FOR CLOSEST POLARIZABLE PAIRS\n');
fprintf('============================================================\n');

if ~isempty(pairs)
    nShow = min(12, size(pairs,1));
    fprintf('   i    j      shortest-image-dist(A)      ||dTij||_F / ||Tij_real||_F\n');

    alphaE = res_thole.Tmeta.alpha;
    rcutE  = res_thole.Tmeta.rcut;

    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice.';
    else
        H = sys.lattice.';
    end

    aH = H(:,1); bH = H(:,2); cH = H(:,3);
    nxmax = ceil(rcutE / norm(aH)) + 1;
    nymax = ceil(rcutE / norm(bH)) + 1;
    nzmax = ceil(rcutE / norm(cH)) + 1;

    for p = 1:nShow
        i = pairs(p,1);
        j = pairs(p,2);

        ri = pos(i,:);
        rj = pos(j,:);

        best_xvec = [];
        best_dist = inf;

        for nx = -nxmax:nxmax
            for ny = -nymax:nymax
                for nz = -nzmax:nzmax
                    if i == j && nx == 0 && ny == 0 && nz == 0
                        continue;
                    end

                    nvec = [nx; ny; nz];
                    Rimg = (H * nvec).';
                    xvec = (ri - rj) + Rimg;
                    dist = norm(xvec);

                    if dist == 0 || dist > rcutE
                        continue;
                    end

                    if dist < best_dist
                        best_dist = dist;
                        best_xvec = xvec;
                    end
                end
            end
        end

        if isempty(best_xvec)
            fprintf('%4d %4d   no image within rcut\n', i, j);
            continue;
        end

        Tbare_pair = ewald.real_space_tensor_block_triclinic(best_xvec, alphaE);
        dT_pair = ewald.real_space_tensor_block_triclinic_thole_correction( ...
            best_xvec, alpha_site(i), alpha_site(j), model.thole_a);

        nb = norm(Tbare_pair, 'fro');
        nd = norm(dT_pair, 'fro');

        if nb > 0
            rel = nd / nb;
        else
            rel = NaN;
        end

        fprintf('%4d %4d   %20.8f      %14.6e\n', i, j, best_dist, rel);
    end
end

fprintf('\n============================================================\n');
fprintf('OPTIONAL STRONGER-THOLE SENSITIVITY CHECK\n');
fprintf('============================================================\n');

model_strong = model;
model_strong.thole_a = 2.0;

params_strong = util.default_params();
params_strong.ewald.mode = 'periodic_triclinic';
params_strong.ewald.auto = true;
params_strong.ewald.tol = 1e-8;
params_strong.ewald.boundary = 'tinfoil';
params_strong.ewald.rcut_fraction = 0.9;
params_strong.ewald.use_thole_real_space = true;
params_strong.ewald.thole_a = model_strong.thole_a;
params_strong.scf.solver = 'direct';

res_strong = calc.run_polarization_calc(crystal, model_strong, opts, params_strong);

fprintf('Strong-Thole total energy = %+ .15e\n', res_strong.energy.total);
fprintf('|delta U| vs bare         = %.3e\n', abs(res_strong.energy.total - res_bare.energy.total));
fprintf('|delta U| vs thole(0.39)  = %.3e\n', abs(res_strong.energy.total - res_thole.energy.total));

dmu_strong = res_strong.mu - res_bare.mu;
fprintf('max |delta mu_i| strong-vs-bare = %.16e\n', max(sqrt(sum(dmu_strong.^2,2))));