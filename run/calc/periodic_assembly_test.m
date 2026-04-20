%% periodic_assembly_toy_regression_test
% Regression test for refactored periodic Ewald assembly on a tiny toy system.
%
% This uses a deliberately simple periodic system:
%   - simple cubic lattice
%   - a few polarizable sites in the unit cell
%   - no chemistry/build pipeline required
%
% The goal is only to verify that:
%   new cache-based periodic assembly == local legacy reference assembly
%
% Notes
%   - This is a correctness test, not a performance benchmark.
%   - The local legacy reference is intentionally implemented inside this file
%     so the minimal refactor branch does not need old source files.

clear; clc; close all;

fprintf('=== periodic toy assembly regression test ===\n');

%% ------------------------------------------------------------------------
% User controls

cfg = struct();

cfg.softening = 0.0;
cfg.use_thole = true;

cfg.ewald_tol = 1e-4;
cfg.boundary = 'tinfoil';
cfg.rcut_fraction = 0.7;

cfg.verbose = true;

% Tiny cubic toy cell
cfg.L_ang = 8.0;

% 4-site toy basis in Angstrom (converted below to bohr)
cfg.site_pos_ang = [ ...
    1.0, 1.0, 1.0; ...
    5.0, 1.0, 1.0; ...
    1.0, 5.0, 1.0; ...
    1.0, 1.0, 5.0];

cfg.site_alpha = [1.20; 1.10; 1.00; 0.90];
cfg.site_is_polarizable = true(4,1);
cfg.thole_a = 0.39;

tol = struct();
tol.Tfro = 1e-10;
tol.Tmax = 1e-10;

%% ------------------------------------------------------------------------
% Build tiny toy periodic polarization system directly in canonical atomic units

ANG2BOHR = 1.8897261254578281;

polsys = struct();
polsys.n_sites = size(cfg.site_pos_ang, 1);
polsys.site_pos = cfg.site_pos_ang * ANG2BOHR;        % bohr
polsys.site_alpha = cfg.site_alpha;                   % atomic units
polsys.site_is_polarizable = cfg.site_is_polarizable;
polsys.thole_a = cfg.thole_a;
polsys.site_charge = zeros(polsys.n_sites, 1);        % elementary charge
polsys.lattice = (cfg.L_ang * ANG2BOHR) * eye(3);     % bohr

polsys.units = struct();
polsys.units.length = 'bohr';
polsys.units.alpha = 'atomic_unit';
polsys.units.charge = 'elementary_charge';

io.assert_atomic_units(polsys);

fprintf('Toy periodic polarization system:\n');
fprintf('  n_sites             = %d\n', polsys.n_sites);
fprintf('  n_polarizable_sites = %d\n', nnz(polsys.site_is_polarizable));
fprintf('  box length          = %.6f bohr\n', cfg.L_ang * ANG2BOHR);

%% ------------------------------------------------------------------------
% External field (zero for pure operator regression)

Eext = zeros(polsys.n_sites, 3);

%% ------------------------------------------------------------------------
% Problem setup and Ewald parameters

scfParams = struct();
scfParams.softening = cfg.softening;
scfParams.use_thole = cfg.use_thole;
scfParams.verbose   = cfg.verbose;
scfParams.printEvery = 25;

problem = thole.prepare_scf_problem(polsys, Eext, scfParams);

H = local_get_direct_lattice(polsys);
ewaldParams = ewald.choose_ewald_params(H, cfg.ewald_tol, cfg.boundary, cfg.rcut_fraction);

fprintf('\nChosen Ewald parameters:\n');
fprintf('  alpha    = %.6g\n', ewaldParams.alpha);
fprintf('  rcut     = %.6g\n', ewaldParams.rcut);
fprintf('  kcut     = %.6g\n', ewaldParams.kcut);
fprintf('  boundary = %s\n', ewaldParams.boundary);

%% ------------------------------------------------------------------------
% New refactored assembler

fprintf('\n--- new refactored periodic assembly ---\n');
tic;
[T_new, parts_new, opinfo_new] = ewald.assemble_periodic_interaction_matrix( ...
    polsys, problem, ewaldParams, scfParams);
time_new = toc;

fprintf('new assembly time = %.6f s\n', time_new);
fprintf('new operator size = %d x %d\n', size(T_new,1), size(T_new,2));
fprintf('new nRealInteractions = %d\n', opinfo_new.nRealInteractions);
fprintf('new num_kvec          = %d\n', opinfo_new.num_kvec);

%% ------------------------------------------------------------------------
% Legacy reference assembler (local implementation)

fprintf('\n--- legacy reference periodic assembly ---\n');
tic;
[T_ref, parts_ref, meta_ref] = local_legacy_periodic_reference( ...
    polsys, problem, ewaldParams, scfParams);
time_ref = toc;

fprintf('reference assembly time = %.6f s\n', time_ref);
fprintf('reference operator size = %d x %d\n', size(T_ref,1), size(T_ref,2));
fprintf('reference num_kvec      = %d\n', meta_ref.num_kvec);

%% ------------------------------------------------------------------------
% Comparisons

D = T_new - T_ref;

Tfro = norm(D, 'fro');
Tmax = max(abs(D(:)));

Dreal  = norm(parts_new.real  - parts_ref.real,  'fro');
Drecip = norm(parts_new.recip - parts_ref.recip, 'fro');
Dself  = norm(parts_new.self  - parts_ref.self,  'fro');
Dsurf  = norm(parts_new.surf  - parts_ref.surf,  'fro');

fprintf('\n--- comparisons ---\n');
fprintf('||T_new - T_ref||_F       = %.16e\n', Tfro);
fprintf('max |T_new - T_ref|       = %.16e\n', Tmax);
fprintf('||real_new  - real_ref||  = %.16e\n', Dreal);
fprintf('||recip_new - recip_ref|| = %.16e\n', Drecip);
fprintf('||self_new  - self_ref||  = %.16e\n', Dself);
fprintf('||surf_new  - surf_ref||  = %.16e\n', Dsurf);

%% ------------------------------------------------------------------------
% Assertions

if Tfro > tol.Tfro
    error('Periodic operator Frobenius mismatch too large: %.3e', Tfro);
end

if Tmax > tol.Tmax
    error('Periodic operator max-entry mismatch too large: %.3e', Tmax);
end

fprintf('\nPASS: refactored periodic assembly matches local legacy reference.\n');

%% ========================================================================
% Local helpers
%% ========================================================================

function [Tpol, parts, meta] = local_legacy_periodic_reference(sys, problem, ewaldParams, scfParams)
%LOCAL_LEGACY_PERIODIC_REFERENCE
% Local reference version based on the old monolithic periodic assembly,
% but restricted to polarizable-only active-space output so it is directly
% comparable to the refactored assembler.

    if nargin < 4 || isempty(scfParams)
        scfParams = struct();
    end

    H = local_get_direct_lattice(sys);

    alpha = ewaldParams.alpha;
    rcut  = ewaldParams.rcut;
    kcut  = ewaldParams.kcut;

    boundary = 'tinfoil';
    if isfield(ewaldParams, 'boundary') && ~isempty(ewaldParams.boundary)
        boundary = ewaldParams.boundary;
    end

    use_thole_real_space = true;
    if isfield(scfParams, 'use_thole') && ~isempty(scfParams.use_thole)
        use_thole_real_space = logical(scfParams.use_thole);
    end

    sites = problem.activeSites(:);
    nPol = problem.nPolSites;
    site_alpha = sys.site_alpha(:);

    Treal  = zeros(3*nPol, 3*nPol);
    Trecip = zeros(3*nPol, 3*nPol);
    Tself  = zeros(3*nPol, 3*nPol);
    Tsurf  = zeros(3*nPol, 3*nPol);

    [kvecs, kmeta] = ewald.enumerate_kvecs_triclinic(H, kcut);

    a = H(:,1);
    b = H(:,2);
    c = H(:,3);

    nxmax = ceil(rcut / norm(a)) + 1;
    nymax = ceil(rcut / norm(b)) + 1;
    nzmax = ceil(rcut / norm(c)) + 1;

    for ai = 1:nPol
        i = sites(ai);
        ii = util.block3(ai);
        ri = sys.site_pos(i, :);

        Tself(ii, ii) = ewald.self_tensor_block_dipole(alpha);

        for aj = 1:nPol
            j = sites(aj);
            jj = util.block3(aj);
            rj = sys.site_pos(j, :);

            Tsurf(ii, jj) = ewald.surface_tensor_block_dipole(H, boundary);

            Trecip(ii, jj) = local_reciprocal_space_tensor_block_triclinic( ...
                ri, rj, H, alpha, kcut, kvecs);

            Tij_real = zeros(3,3);
            rij0 = ri - rj;

            for nx = -nxmax:nxmax
                for ny = -nymax:nymax
                    for nz = -nzmax:nzmax
                        if i == j && nx == 0 && ny == 0 && nz == 0
                            continue;
                        end

                        nvec = [nx; ny; nz];
                        Rimg = (H * nvec).';
                        xvec = rij0 + Rimg;
                        x = norm(xvec);

                        if x == 0 || x > rcut
                            continue;
                        end

                        Tij_block = local_real_space_tensor_block_triclinic(xvec, alpha);

                        if use_thole_real_space
                            Tij_block = Tij_block + ...
                                local_real_space_tensor_block_triclinic_thole_correction( ...
                                    xvec, site_alpha(i), site_alpha(j), sys.thole_a);
                        end

                        Tij_real = Tij_real + Tij_block;
                    end
                end
            end

            Treal(ii, jj) = Tij_real;
        end
    end

    Tpol = Treal + Trecip + Tself + Tsurf;

    parts = struct();
    parts.real  = Treal;
    parts.recip = Trecip;
    parts.self  = Tself;
    parts.surf  = Tsurf;

    meta = struct();
    meta.real_image_bounds = [nxmax nymax nzmax];
    meta.num_kvec = kmeta.num_kvec;
    meta.volume = abs(det(H));
    meta.boundary = boundary;
    meta.kcut = kcut;
    meta.rcut = rcut;
    meta.alpha = alpha;
end

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

function Tij = local_reciprocal_space_tensor_block_triclinic(ri, rj, H, alpha, kcut, kvecs)
%LOCAL_RECIPROCAL_SPACE_TENSOR_BLOCK_TRICLINIC
% Legacy-style reciprocal-space Ewald dipole tensor block.

    if nargin < 6 || isempty(kvecs)
        [kvecs, ~] = ewald.enumerate_kvecs_triclinic(H, kcut);
    end

    V = abs(det(H));
    rij = ri - rj;
    Tij = zeros(3, 3);

    for m = 1:size(kvecs, 1)
        kvec = kvecs(m, :);
        k2 = dot(kvec, kvec);

        pref = (2 * pi / V) * exp(-k2 / (4 * alpha^2)) / k2;
        phase = dot(kvec, rij);

        kkT = kvec(:) * kvec(:).';
        Tij = Tij + 2 * pref * cos(phase) * kkT;
    end
end

function Tij = local_real_space_tensor_block_triclinic(xvec, alpha)
%LOCAL_REAL_SPACE_TENSOR_BLOCK_TRICLINIC
% Legacy-style real-space Ewald dipole tensor block.

    x = norm(xvec);

    if x == 0
        Tij = zeros(3,3);
        return;
    end

    erfcax = erfc(alpha * x);
    expax2 = exp(-(alpha^2) * (x^2));

    B = erfcax / x^3 + (2 * alpha / sqrt(pi)) * expax2 / x^2;
    C = 3 * erfcax / x^5 + (2 * alpha / sqrt(pi)) * ...
        (2 * alpha^2 / x^2 + 3 / x^4) * expax2;

    xxT = xvec(:) * xvec(:).';
    Tij = B * eye(3) - C * xxT;
end

function dTij = local_real_space_tensor_block_triclinic_thole_correction(xvec, alpha_i, alpha_j, thole_a)
%LOCAL_REAL_SPACE_TENSOR_BLOCK_TRICLINIC_THOLE_CORRECTION
% Legacy-style Thole correction added to the bare Ewald real-space tensor.

    x = norm(xvec);

    if x == 0
        dTij = zeros(3,3);
        return;
    end

    tf = thole.thole_f3f5_factors(x, alpha_i, alpha_j, thole_a);

    xxT = xvec(:) * xvec(:).';
    dTij = (tf.l3 / x^3) * eye(3) - (3 * tf.l5 / x^5) * xxT;
end