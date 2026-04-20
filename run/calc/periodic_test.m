%% periodic_vs_nonperiodic_operator_limit_test
% Compare assembled periodic and nonperiodic operators for a fixed finite
% fragment as vacuum padding is increased.
%
% Goal:
%   Determine whether the periodic operator itself approaches the
%   nonperiodic operator in the large-vacuum limit.
%
% This is a more direct diagnostic than comparing solved dipoles alone.

clear; clc; close all;

fprintf('=== periodic vs nonperiodic operator vacuum-limit test ===\n');

%% ------------------------------------------------------------------------
% User controls

cfg = struct();

cfg.softening = 0.0;
cfg.use_thole = true;
cfg.boundary = 'tinfoil';

cfg.nonperiodic_softening = 0.0;
cfg.periodic_ewald_tol = 1e-6;
cfg.periodic_rcut_fraction = 0.7;

cfg.verbose = true;

% External field is not used for operator comparison, but prepare_scf_problem
% wants one. Keep it zero here.
cfg.Eext = [0.0, 0.0, 0.0];

% Vacuum padding sweep in Angstrom
cfg.box_list_ang = [10, 100, 200, 300, 400, 500];

% Small fixed finite fragment in Angstrom, converted below to bohr
cfg.site_pos_ang = [ ...
    0.0,  0.0,  0.0; ...
    3.1,  0.2, -0.1; ...
    0.5,  3.0,  0.3; ...
   -0.3,  0.6,  2.9];

cfg.site_alpha = [1.20; 1.10; 1.00; 0.90];
cfg.site_is_polarizable = true(4,1);
cfg.thole_a = 0.39;

%% ------------------------------------------------------------------------
% Build canonical finite fragment in atomic units

ANG2BOHR = 1.8897261254578281;

sys_np = struct();
sys_np.n_sites = size(cfg.site_pos_ang, 1);
sys_np.site_pos = cfg.site_pos_ang * ANG2BOHR;          % bohr
sys_np.site_alpha = cfg.site_alpha;                     % atomic units
sys_np.site_is_polarizable = cfg.site_is_polarizable;
sys_np.thole_a = cfg.thole_a;
sys_np.site_charge = zeros(sys_np.n_sites, 1);          % elementary charge

sys_np.units = struct();
sys_np.units.length = 'bohr';
sys_np.units.alpha = 'atomic_unit';
sys_np.units.charge = 'elementary_charge';

io.assert_atomic_units(sys_np);

fprintf('Finite fragment:\n');
fprintf('  n_sites             = %d\n', sys_np.n_sites);
fprintf('  n_polarizable_sites = %d\n', nnz(sys_np.site_is_polarizable));

%% ------------------------------------------------------------------------
% Dummy external field for problem preparation

Eext = repmat(cfg.Eext, sys_np.n_sites, 1);

%% ------------------------------------------------------------------------
% Nonperiodic assembled reference operator

fprintf('\n--- nonperiodic assembled reference ---\n');

scfNP = struct();
scfNP.softening = cfg.nonperiodic_softening;
scfNP.use_thole = cfg.use_thole;
scfNP.verbose = cfg.verbose;
scfNP.printEvery = 25;

problem_np = thole.prepare_scf_problem(sys_np, Eext, scfNP);

tic;
Tnp = ewald.assemble_nonperiodic_interaction_matrix( ...
    sys_np, problem_np, struct(), scfNP);
time_np_asm = toc;

fprintf('nonperiodic assembly time = %.6f s\n', time_np_asm);
fprintf('nonperiodic operator size = %d x %d\n', size(Tnp,1), size(Tnp,2));

fprintf('\nNonperiodic diagonal blocks:\n');
for a = 1:problem_np.nPolSites
    Ia = util.block3(a);
    fprintf('site %d:\n', a);
    disp(Tnp(Ia, Ia));
end

%% ------------------------------------------------------------------------
% Periodic vacuum-padding sweep

fprintf('\n--- periodic operator vacuum-padding sweep ---\n');

nBox = numel(cfg.box_list_ang);

results = struct();
results.L_ang = cfg.box_list_ang(:);
results.L_bohr = cfg.box_list_ang(:) * ANG2BOHR;
results.Tfro = nan(nBox, 1);
results.Tmax = nan(nBox, 1);
results.diag_fro = nan(nBox, 1);
results.offdiag_fro = nan(nBox, 1);
results.assembly_time = nan(nBox, 1);
results.nRealInteractions = nan(nBox, 1);
results.num_kvec = nan(nBox, 1);

for k = 1:nBox
    L_ang = cfg.box_list_ang(k);
    L_bohr = L_ang * ANG2BOHR;

    fprintf('\n>>> cubic box length = %.3f A (%.6f bohr)\n', L_ang, L_bohr);

    sys_per = sys_np;
    sys_per.lattice = L_bohr * eye(3);

    scfPer = struct();
    scfPer.softening = cfg.softening;
    scfPer.use_thole = cfg.use_thole;
    scfPer.verbose = cfg.verbose;
    scfPer.printEvery = 25;

    problem_per = thole.prepare_scf_problem(sys_per, Eext, scfPer);

    H = local_get_direct_lattice(sys_per);
    ewaldParams = ewald.choose_ewald_params(H, ...
        cfg.periodic_ewald_tol, cfg.boundary, cfg.periodic_rcut_fraction);

    fprintf('  alpha = %.6g | rcut = %.6g | kcut = %.6g\n', ...
        ewaldParams.alpha, ewaldParams.rcut, ewaldParams.kcut);

    tAsm = tic;
    [Tper, ~, opinfo_per] = ewald.assemble_periodic_interaction_matrix( ...
        sys_per, problem_per, ewaldParams, scfPer);
    results.assembly_time(k) = toc(tAsm);

    D = Tper - Tnp;

    results.Tfro(k) = norm(D, 'fro');
    results.Tmax(k) = max(abs(D(:)));
    results.diag_fro(k) = local_diag_block_fro(D, problem_np.nPolSites);
    results.offdiag_fro(k) = local_offdiag_block_fro(D, problem_np.nPolSites);
    results.nRealInteractions(k) = opinfo_per.nRealInteractions;
    results.num_kvec(k) = opinfo_per.num_kvec;

    fprintf('  ||Tper - Tnp||_F    = %.16e\n', results.Tfro(k));
    fprintf('  max |Tper - Tnp|    = %.16e\n', results.Tmax(k));
    fprintf('  diag-block Fro norm = %.16e\n', results.diag_fro(k));
    fprintf('  offdiag Fro norm    = %.16e\n', results.offdiag_fro(k));
    fprintf('  nRealInteractions   = %d\n', opinfo_per.nRealInteractions);
    fprintf('  num_kvec            = %d\n', opinfo_per.num_kvec);

    if k == 1 || k == nBox
        fprintf('\n  Block comparison for L = %.3f A\n', L_ang);
        for a = 1:problem_np.nPolSites
            Ia = util.block3(a);
            fprintf('  diagonal block site %d: ||D_ii||_F = %.16e\n', ...
                a, norm(D(Ia, Ia), 'fro'));
        end

        if problem_np.nPolSites >= 2
            I1 = util.block3(1);
            I2 = util.block3(2);
            fprintf('  sample offdiag block (1,2): ||D_12||_F = %.16e\n', ...
                norm(D(I1, I2), 'fro'));
            fprintf('  sample offdiag block difference D_12 = \n');
            disp(D(I1, I2));
        end
    end
end

%% ------------------------------------------------------------------------
% Summary

fprintf('\n--- operator vacuum-limit summary ---\n');
fprintf('%10s  %14s  %14s  %14s  %14s  %12s  %12s\n', ...
    'L(A)', '||D||_F', 'max|D|', 'diag_F', 'offdiag_F', 'nReal', 'nK');

for k = 1:nBox
    fprintf('%10.3f  %14.6e  %14.6e  %14.6e  %14.6e  %12d  %12d\n', ...
        results.L_ang(k), results.Tfro(k), results.Tmax(k), ...
        results.diag_fro(k), results.offdiag_fro(k), ...
        results.nRealInteractions(k), results.num_kvec(k));
end

fprintf('\nInterpretation:\n');
fprintf('  If periodic -> nonperiodic in the large-vacuum limit, the operator\n');
fprintf('  difference metrics above should trend toward zero.\n');
fprintf('Done.\n');

%% ------------------------------------------------------------------------
% Targeted off-diagonal sign check at the largest box

fprintf('\n--- targeted off-diagonal sign check (largest box) ---\n');

L_ang = cfg.box_list_ang(end);
L_bohr = L_ang * ANG2BOHR;

sys_per = sys_np;
sys_per.lattice = L_bohr * eye(3);

scfPer = struct();
scfPer.softening = cfg.softening;
scfPer.use_thole = cfg.use_thole;
scfPer.verbose = false;
scfPer.printEvery = 25;

problem_per = thole.prepare_scf_problem(sys_per, Eext, scfPer);

H = local_get_direct_lattice(sys_per);
ewaldParams = ewald.choose_ewald_params(H, ...
    cfg.periodic_ewald_tol, cfg.boundary, cfg.periodic_rcut_fraction);

[Tper, ~, ~] = ewald.assemble_periodic_interaction_matrix( ...
    sys_per, problem_per, ewaldParams, scfPer);

I1 = util.block3(1);
I2 = util.block3(2);

Tnp12 = Tnp(I1, I2);
Tper12 = Tper(I1, I2);

fprintf('Largest box L = %.3f A\n', L_ang);

fprintf('\nTnp(1,2) =\n');
disp(Tnp12);

fprintf('Tper(1,2) =\n');
disp(Tper12);

fprintf('Tper(1,2) - Tnp(1,2) =\n');
disp(Tper12 - Tnp12);

fprintf('Tper(1,2) + Tnp(1,2) =\n');
disp(Tper12 + Tnp12);

fprintf('||Tper12 - Tnp12||_F = %.16e\n', norm(Tper12 - Tnp12, 'fro'));
fprintf('||Tper12 + Tnp12||_F = %.16e\n', norm(Tper12 + Tnp12, 'fro'));

%% ========================================================================
% Local helpers
%% ========================================================================

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('local_get_direct_lattice:MissingLattice', ...
            'Need sys.super_lattice or sys.lattice.');
    end
end

function v = local_diag_block_fro(M, nPol)
    acc = 0.0;
    for a = 1:nPol
        Ia = util.block3(a);
        acc = acc + norm(M(Ia, Ia), 'fro')^2;
    end
    v = sqrt(acc);
end

function v = local_offdiag_block_fro(M, nPol)
    acc = 0.0;
    for a = 1:nPol
        Ia = util.block3(a);
        for b = 1:nPol
            if a == b
                continue;
            end
            Ib = util.block3(b);
            acc = acc + norm(M(Ia, Ib), 'fro')^2;
        end
    end
    v = sqrt(acc);
end