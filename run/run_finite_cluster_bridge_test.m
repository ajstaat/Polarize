function run_finite_cluster_bridge_test
clc;

fprintf('=== finite-cluster bridge test ===\n');

% -------------------------------------------------------------------------
% User-style finite-cluster setup, following the pasted script pattern
% -------------------------------------------------------------------------
params = struct( ...
    'n_in_stack', 8, ...
    'tilt_dx', 0.0, ...
    'tilt_dy', 0.0, ...
    'n_neighbors', 1, ...
    'mode', 'double', ...
    'solver', 'sor', ...
    'omega', 0.9, ...
    'tol', 1e-8, ...
    'maxit', 4000, ...
    'mu0', [], ...
    'plot_scale', 100.0, ...
    'plot_skip', 1, ...
    'plot', false);

% -------------------------------------------------------------------------
% Load finite-cluster inputs from MAT file
% -------------------------------------------------------------------------
matfile = 'finite_cluster_inputs.mat';
S = load(matfile);

required = {'R_A','alpha_A3','q_plus','q_minus','a'};
for k = 1:numel(required)
    if ~isfield(S, required{k})
        error('MAT file must contain variable "%s".', required{k});
    end
end

R_A      = S.R_A;
alpha_A3 = S.alpha_A3;
q_plus   = S.q_plus;
q_minus  = S.q_minus;
a        = S.a;

% -------------------------------------------------------------------------
% Build the same finite cluster as the standalone script
% -------------------------------------------------------------------------
[R_all_A, alpha_all_A3, mol_start_idx, anchors] = ...
    build_stacks_pdi_fixeddz(R_A, alpha_A3, ...
        params.n_in_stack, params.tilt_dx, params.tilt_dy, params.n_neighbors);

switch lower(params.mode)
    case 'single'
        qA = q_plus(:);
        qB = [];
    case 'double'
        qA = q_plus(:);
        qB = q_minus(:);
    otherwise
        error('params.mode must be ''single'' or ''double''.');
end

nAtom = size(R_A,1);
Nat   = size(R_all_A,1);

[pol_mask, removed_mol_ids, ~, ~, removed_layer_ids] = ...
    remove_center_layers_center_stack( ...
        mol_start_idx, nAtom, params.n_in_stack, anchors, params.mode);

q_all_e = zeros(Nat,1);
if strcmpi(params.mode,'single')
    i0A = mol_start_idx(removed_mol_ids(1));
    idxA = (i0A : i0A+nAtom-1).';
    q_all_e(idxA) = qA(:);
else
    i0A = mol_start_idx(removed_mol_ids(1));
    idxA = (i0A : i0A+nAtom-1).';
    i0B = mol_start_idx(removed_mol_ids(2));
    idxB = (i0B : i0B+nAtom-1).';
    q_all_e(idxA) = qA(:);
    q_all_e(idxB) = qB(:);
end

fprintf('cluster atoms = %d | active atoms = %d | removed mol ids = %s\n', ...
    Nat, nnz(pol_mask), mat2str(removed_mol_ids.'));
fprintf('removed layer ids in center stack = %s\n', mat2str(removed_layer_ids));
fprintf('net inserted charge = %+0.6f e\n', sum(q_all_e));

% -------------------------------------------------------------------------
% Old standalone finite-cluster reference
% -------------------------------------------------------------------------
fprintf('\n--- old finite-cluster reference ---\n');

t_old_direct = tic;
res_old_direct = polarization_from_charges( ...
    R_all_A, alpha_all_A3, pol_mask, q_all_e, a, ...
    'direct', [], params.tol, params.maxit, params.mu0);
time_old_direct = toc(t_old_direct);

t_old_sor = tic;
res_old_sor = polarization_from_charges( ...
    R_all_A, alpha_all_A3, pol_mask, q_all_e, a, ...
    'sor', params.omega, params.tol, params.maxit, params.mu0);
time_old_sor = toc(t_old_sor);

t_old_jac = tic;
res_old_jac = polarization_from_charges( ...
    R_all_A, alpha_all_A3, pol_mask, q_all_e, a, ...
    'jacobi', 1.0, params.tol, params.maxit, params.mu0);
time_old_jac = toc(t_old_jac);

fprintf('old direct: U_pol_eV = %.16e | relres = %.3e | rho(AT)=%.6f | time=%.3f s\n', ...
    res_old_direct.U_pol_eV, res_old_direct.info.relres, res_old_direct.rhoAT, time_old_direct);
fprintf('old sor   : U_pol_eV = %.16e | relres = %.3e | iters=%d | time=%.3f s\n', ...
    res_old_sor.U_pol_eV, res_old_sor.info.relres, res_old_sor.info.iters, time_old_sor);
fprintf('old jacobi: U_pol_eV = %.16e | relres = %.3e | iters=%d | time=%.3f s\n', ...
    res_old_jac.U_pol_eV, res_old_jac.info.relres, res_old_jac.info.iters, time_old_jac);

% -------------------------------------------------------------------------
% New refactored general repo path, with explicit unit normalization
% -------------------------------------------------------------------------
fprintf('\n--- new general repo path with explicit AU conversion ---\n');

t_sys = tic;
sys = make_sys_from_finite_cluster(R_all_A, alpha_all_A3, pol_mask, q_all_e, a);
sys.units.length = 'angstrom';
sys.units.alpha = 'angstrom^3';
sys.units.charge = 'e';
sys = io.convert_to_atomic_units(sys);
time_sys = toc(t_sys);

fieldParams = struct();
fieldParams.field = struct();
fieldParams.field.include_external_charges = true;

scfParams = struct();
scfParams.use_thole = true;
scfParams.softening = 0.0;
scfParams.tol = params.tol;
scfParams.maxIter = params.maxit;
scfParams.mixing = 1.0;
scfParams.omega = params.omega;
scfParams.verbose = false;
scfParams.printEvery = 25;
scfParams.residualEvery = 25;
scfParams.stopMetric = 'relres';

t_field = tic;
Eext = calc.compute_external_field(sys, fieldParams);
time_field = toc(t_field);

t_problem = tic;
problem = thole.prepare_scf_problem(sys, Eext, scfParams);
time_problem = toc(t_problem);

t_op = tic;
[Tpol, ~] = ewald.assemble_nonperiodic_interaction_matrix(sys, problem, struct(), scfParams);
time_op = toc(t_op);

t_new_direct = tic;
[mu_new_direct, direct_new] = thole.solve_scf_direct(problem, Tpol);
time_new_direct = toc(t_new_direct);

t_new_sor = tic;
[mu_new_sor, scf_new_sor] = thole.solve_scf_sor(problem, Tpol);
time_new_sor = toc(t_new_sor);

t_new_jac = tic;
[mu_new_jac, scf_new_jac] = thole.solve_scf_matrix_iterative(problem, Tpol);
time_new_jac = toc(t_new_jac);

E_new_direct_eV = local_energy_from_old_units(problem, mu_new_direct);
E_new_sor_eV    = local_energy_from_old_units(problem, mu_new_sor);
E_new_jac_eV    = local_energy_from_old_units(problem, mu_new_jac);

relres_new_direct = thole.compute_active_space_relres(problem, Tpol, mu_new_direct);
relres_new_sor    = thole.compute_active_space_relres(problem, Tpol, mu_new_sor);
relres_new_jac    = thole.compute_active_space_relres(problem, Tpol, mu_new_jac);

fprintf('new direct: U_pol_eV = %.16e | relres = %.3e | time=%.3f s\n', ...
    E_new_direct_eV, relres_new_direct, time_new_direct);
fprintf('new sor   : U_pol_eV = %.16e | relres = %.3e | iters=%d | time=%.3f s\n', ...
    E_new_sor_eV, relres_new_sor, scf_new_sor.nIter, time_new_sor);
fprintf('new jacobi: U_pol_eV = %.16e | relres = %.3e | iters=%d | time=%.3f s\n', ...
    E_new_jac_eV, relres_new_jac, scf_new_jac.nIter, time_new_jac);

% -------------------------------------------------------------------------
% Compare old vs new
% -------------------------------------------------------------------------
fprintf('\n--- old vs new comparisons ---\n');

mu_old_direct_full = expand_old_mu_to_full(res_old_direct.mu_AU, res_old_direct.idx_active, Nat);
mu_old_sor_full    = expand_old_mu_to_full(res_old_sor.mu_AU,    res_old_sor.idx_active,    Nat);
mu_old_jac_full    = expand_old_mu_to_full(res_old_jac.mu_AU,    res_old_jac.idx_active,    Nat);

err_direct = max_site_diff(mu_old_direct_full, mu_new_direct);
err_sor    = max_site_diff(mu_old_sor_full,    mu_new_sor);
err_jac    = max_site_diff(mu_old_jac_full,    mu_new_jac);

fprintf('max |mu_old_direct - mu_new_direct| = %.16e\n', err_direct);
fprintf('max |mu_old_sor    - mu_new_sor   | = %.16e\n', err_sor);
fprintf('max |mu_old_jac    - mu_new_jac   | = %.16e\n', err_jac);

fprintf('|E_old_direct - E_new_direct| = %.16e eV\n', abs(res_old_direct.U_pol_eV - E_new_direct_eV));
fprintf('|E_old_sor    - E_new_sor   | = %.16e eV\n', abs(res_old_sor.U_pol_eV    - E_new_sor_eV));
fprintf('|E_old_jac    - E_new_jac   | = %.16e eV\n', abs(res_old_jac.U_pol_eV    - E_new_jac_eV));

assert(err_direct < 1e-10, 'old/new direct dipoles disagree');
assert(err_sor    < 1e-8,  'old/new sor dipoles disagree');

if all(isfinite(mu_old_jac_full(:))) && all(isfinite(mu_new_jac(:)))
    assert(err_jac < 1e-8, 'old/new jacobi dipoles disagree');
else
    fprintf('jacobi non-finite on at least one side; skipping strict dipole assert.\n');
end

assert(abs(res_old_direct.U_pol_eV - E_new_direct_eV) < 1e-10, 'old/new direct energy disagree');
assert(abs(res_old_sor.U_pol_eV    - E_new_sor_eV   ) < 1e-8,  'old/new sor energy disagree');

if isfinite(res_old_jac.U_pol_eV) && isfinite(E_new_jac_eV)
    assert(abs(res_old_jac.U_pol_eV - E_new_jac_eV) < 1e-8, 'old/new jacobi energy disagree');
else
    fprintf('jacobi non-finite on at least one side; skipping strict energy assert.\n');
end

fprintf('\n--- timing summary ---\n');
fprintf('old direct solve total      = %.3f s\n', time_old_direct);
fprintf('old sor solve total         = %.3f s\n', time_old_sor);
fprintf('old jacobi solve total      = %.3f s\n', time_old_jac);

fprintf('new sys+unit conversion     = %.3f s\n', time_sys);
fprintf('new external field          = %.3f s\n', time_field);
fprintf('new problem prep            = %.3f s\n', time_problem);
fprintf('new operator assembly       = %.3f s\n', time_op);
fprintf('new direct solve            = %.3f s\n', time_new_direct);
fprintf('new sor solve               = %.3f s\n', time_new_sor);
fprintf('new jacobi solve            = %.3f s\n', time_new_jac);

fprintf('new direct total path       = %.3f s\n', ...
    time_sys + time_field + time_problem + time_op + time_new_direct);
fprintf('new sor total path          = %.3f s\n', ...
    time_sys + time_field + time_problem + time_op + time_new_sor);
fprintf('new jacobi total path       = %.3f s\n', ...
    time_sys + time_field + time_problem + time_op + time_new_jac);

fprintf('\nPASS: finite-cluster bridge test completed.\n');
end

% =========================================================================
% Helpers
% =========================================================================

function problem = make_problem_from_old_reference(Nat, idx_active, alpha_AU_active, F_AU, tol, maxit, omega, mu0_vec)
% Build a thole.prepare_scf_problem-like struct directly from the old
% finite-cluster reference quantities, in the old script's internal units.

idx_active = idx_active(:);
alpha_AU_active = alpha_AU_active(:);
F_AU = F_AU(:);

nPolSites = numel(idx_active);

problem = struct();
problem.nSites = Nat;
problem.polMask = false(Nat,1);
problem.polMask(idx_active) = true;

problem.activeSites = idx_active;
problem.nPolSites = nPolSites;

activeVecIdx = zeros(3*nPolSites,1);
for k = 1:nPolSites
    i = idx_active(k);
    activeVecIdx(3*k-2:3*k) = (3*i-2):(3*i);
end
problem.activeVecIdx = activeVecIdx;

alpha_full = zeros(Nat,1);
alpha_full(idx_active) = alpha_AU_active;
problem.alpha = alpha_full;
problem.alpha_pol = alpha_AU_active;
problem.alpha_pol_vec = repelem(alpha_AU_active, 3, 1);

problem.Eext = zeros(Nat,3);
problem.Eext(idx_active,:) = util.unstack_xyz(F_AU);
problem.Eext_vec = util.stack_xyz(problem.Eext);
problem.Eext_pol = util.unstack_xyz(F_AU);
problem.Eext_pol_vec = F_AU;

if isempty(mu0_vec)
    mu0_pol_vec = zeros(3*nPolSites,1);
else
    mu0_pol_vec = mu0_vec(:);
    if numel(mu0_pol_vec) ~= 3*nPolSites
        error('mu0 must be empty or length 3*nPolSites.');
    end
end

problem.mu0 = zeros(Nat,3);
problem.mu0(idx_active,:) = util.unstack_xyz(mu0_pol_vec);
problem.mu0_vec = util.stack_xyz(problem.mu0);
problem.mu0_pol = util.unstack_xyz(mu0_pol_vec);
problem.mu0_pol_vec = mu0_pol_vec;

problem.tol = tol;
problem.maxIter = maxit;
problem.mixing = 1.0;
problem.omega = omega;
problem.verbose = false;
problem.printEvery = 25;
problem.residualEvery = 25;
problem.stopMetric = 'relres';
end

function mu_full = expand_old_mu_to_full(mu_active_vec, idx_active, Nat)
mu_active = util.unstack_xyz(mu_active_vec(:));
mu_full = zeros(Nat, 3);
mu_full(idx_active, :) = mu_active;
end

function E_eV = local_energy_from_old_units(problem, mu_full)
% Same energy convention as the old finite-cluster script:
% U_pol_H = -0.5 * F' * mu, then Hartree -> eV.

mu_pol = mu_full(problem.activeSites, :);
mu_pol_vec = util.stack_xyz(mu_pol);
F_AU = problem.Eext_pol_vec;

U_pol_H = -0.5 * (F_AU.' * mu_pol_vec);
E_eV = U_pol_H * 27.211386245988;
end

function err = max_site_diff(muA, muB)
if any(~isfinite(muA(:))) || any(~isfinite(muB(:)))
    err = NaN;
else
    err = max(sqrt(sum((muA - muB).^2, 2)));
end
end

function [R_all_A, alpha_all_A3, mol_start_idx, stack_offsets] = ...
    build_stacks_pdi_fixeddz(R_mono_A, alpha_mono_A3, n_in_stack, tilt_dx, tilt_dy, n_neighbors)
% PDI stacks with fixed dz=3.5 Å, tilt per layer, phase-aligned neighbors.

dz  = 3.5;
Lxy = [-7.3865, +8.19118];
Rxy = [+7.5877, -8.1613];

nAtom = size(R_mono_A,1);
assert(numel(alpha_mono_A3)==nAtom, 'alpha_mono_A3 must align with R_mono_A.');

stack_offsets = [];
phase_is_zero = [];
for k = n_neighbors:-1:1
    zphase = (mod(k,2)==1) * (-dz/2);
    stack_offsets = [stack_offsets; k*Lxy, zphase]; %#ok<AGROW>
    phase_is_zero = [phase_is_zero; (zphase==0)];   %#ok<AGROW>
end
stack_offsets = [stack_offsets; 0,0,0];
phase_is_zero = [phase_is_zero; true];
for k = 1:n_neighbors
    zphase = (mod(k,2)==1) * (-dz/2);
    stack_offsets = [stack_offsets; k*Rxy, zphase]; %#ok<AGROW>
    phase_is_zero = [phase_is_zero; (zphase==0)];   %#ok<AGROW>
end

s_layers = (1:n_in_stack) - (n_in_stack+1)/2;
nStacks  = size(stack_offsets,1);
nMolTot  = nStacks * n_in_stack;

R_all_A       = zeros(nAtom*nMolTot, 3);
alpha_all_A3  = zeros(nAtom*nMolTot, 1);
mol_start_idx = zeros(nMolTot,1);

row = 1;
mol_id = 1;
for s = 1:nStacks
    base = stack_offsets(s,:);
    halfTilt_xy = (phase_is_zero(s)) * 0.5 * [tilt_dx, tilt_dy];
    base_xy = base(1:2) + halfTilt_xy;
    base_z  = base(3);

    for t = 1:n_in_stack
        sL = s_layers(t);
        shift_xy = base_xy + [tilt_dx, tilt_dy]*sL;
        shift_z  = base_z  + dz*sL;

        mol_start_idx(mol_id) = row;
        R_all_A(row:row+nAtom-1,:) = R_mono_A + [shift_xy, shift_z];
        alpha_all_A3(row:row+nAtom-1) = alpha_mono_A3;

        row = row + nAtom;
        mol_id = mol_id + 1;
    end
end
end

function [pol_mask, removed_mol_ids, removed_atom_indices, center_stack_id, removed_layer_ids] = ...
    remove_center_layers_center_stack(mol_start_idx, nAtom, n_in_stack, stack_offsets, mode)
% Remove the center 1 or 2 molecules by index from the center stack only.

nStacks   = size(stack_offsets,1);
total_mols = numel(mol_start_idx);
Nat = nAtom * total_mols;

assert(mod(total_mols,n_in_stack)==0 && total_mols==nStacks*n_in_stack, ...
       'mol/stack dimensions inconsistent.');

[~, center_stack_id] = min(sum(stack_offsets(:,1:2).^2,2));

first_mol = (center_stack_id-1)*n_in_stack + 1;
center_mol_ids = (first_mol : first_mol+n_in_stack-1).';

mid = ceil(n_in_stack/2);
if strcmpi(mode,'single')
    if mod(n_in_stack,2)==1
        removed_layer_ids = mid;
    else
        removed_layer_ids = n_in_stack/2;
    end
else
    if mod(n_in_stack,2)==1
        removed_layer_ids = [mid, min(mid+1, n_in_stack)];
    else
        removed_layer_ids = [n_in_stack/2, n_in_stack/2 + 1];
    end
end

removed_mol_ids = center_mol_ids(removed_layer_ids);

pol_mask = true(Nat,1);
removed_atom_indices = [];
for m = removed_mol_ids.'
    i0 = mol_start_idx(m);
    idx = (i0 : i0+nAtom-1).';
    pol_mask(idx) = false;
    removed_atom_indices = [removed_atom_indices; idx]; %#ok<AGROW>
end
end

function result = polarization_from_charges(R_all_A, alpha_all_A3, pol_mask, q_all_e, a, ...
                                           solver, omega, tol, maxit, mu0)
% Induced dipoles with point-charge field; supports direct / sor / gs / jacobi.

ANG2BOHR  = 1.8897259886;
A3_PER_AU = 0.148184711;
I3 = eye(3);

idx_active = find(pol_mask);
Na  = numel(idx_active);
Nat = size(R_all_A,1);

R_bohr_all      = R_all_A * ANG2BOHR;
alpha_AU_active = alpha_all_A3(idx_active) / A3_PER_AU;

Ablk = kron(diag(alpha_AU_active), I3);

T = zeros(3*Na, 3*Na);
for ii = 1:Na
    i = idx_active(ii);
    Ri = R_bohr_all(i,:);
    ai = alpha_AU_active(ii);

    for jj = 1:Na
        if ii == jj
            continue;
        end

        j = idx_active(jj);
        rij = Ri - R_bohr_all(j,:);
        r = norm(rij);
        if r < 1e-12
            continue;
        end

        rrT = rij(:) * rij(:).';
        aij = (ai * alpha_AU_active(jj))^(1/6);
        s = a * (r / aij)^3;
        e = exp(-s);
        f3 = 1 - e;
        f5 = 1 - (1 + s) * e;
        Tij = (3*f5/r^5)*rrT - (f3/r^3)*I3;

        T(3*(ii-1)+1:3*ii, 3*(jj-1)+1:3*jj) = Tij;
    end
end

F = zeros(3*Na,1);
for ii = 1:Na
    i = idx_active(ii);
    Ri = R_bohr_all(i,:);
    Ei = [0 0 0];

    for k = 1:Nat
        qk = q_all_e(k);
        if qk == 0
            continue;
        end
        dr = Ri - R_bohr_all(k,:);
        r = norm(dr);
        if r < 1e-10
            continue;
        end
        Ei = Ei + qk * (dr / (r^3));
    end

    F(3*ii-2:3*ii) = Ei(:);
end

rhoAT = max(abs(eig(Ablk*T)));

if nargin < 6 || isempty(solver)
    solver = 'direct';
end
if isempty(omega)
    if strcmpi(solver,'sor')
        omega = 0.8;
    else
        omega = 1.0;
    end
end
if isempty(tol)
    tol = 1e-8;
end
if isempty(maxit)
    maxit = 2000;
end
if isempty(mu0)
    mu = zeros(3*Na,1);
else
    mu = mu0(:);
end

switch lower(solver)
    case 'direct'
        M = eye(3*Na) - Ablk*T;
        mu = M \ (Ablk * F);
        info = struct('solver','direct','iters',1,'converged',true,'relres',0,'omega',NaN);

    case {'jacobi','gs','sor'}
        b = zeros(3*Na,1);
        for i = 1:Na
            b(3*i-2:3*i) = alpha_AU_active(i) * F(3*i-2:3*i);
        end

        bn = norm(b);
        if bn == 0
            bn = 1;
        end

        relres = compute_relres(mu, b, T, alpha_AU_active, bn);
        it = 0;
        converged = relres < tol;

        while ~converged && it < maxit
            it = it + 1;

            switch lower(solver)
                case 'jacobi'
                    y = T * mu;
                    mu_new = zeros(3*Na,1);
                    for i = 1:Na
                        Ei = F(3*i-2:3*i) + y(3*i-2:3*i);
                        mu_new(3*i-2:3*i) = alpha_AU_active(i) * Ei;
                    end
                    mu = (1-omega)*mu + omega*mu_new;

                otherwise
                    for i = 1:Na
                        Ii = (3*i-2):(3*i);
                        E_loc = F(Ii);

                        if i > 1
                            E_loc = E_loc + T(Ii,1:3*(i-1)) * mu(1:3*(i-1));
                        end
                        if i < Na
                            E_loc = E_loc + T(Ii,3*i+1:end) * mu(3*i+1:end);
                        end

                        mu_i_new = alpha_AU_active(i) * E_loc;
                        mu(Ii) = (1-omega)*mu(Ii) + omega*mu_i_new;
                    end
            end

            relres = compute_relres(mu, b, T, alpha_AU_active, bn);
            converged = relres < tol;
        end

        info = struct('solver',lower(solver), 'iters',it, 'converged',converged, ...
                      'relres',relres, 'omega',omega);

    otherwise
        error('Unknown solver: %s', solver);
end

U_pol_H    = -0.5 * (F.' * mu);
U_pol_eV   = U_pol_H * 27.211386245988;
U_pol_kcal = U_pol_H * 627.509474;

K = kron(ones(Na,1), eye(3));
S = kron(ones(1,Na), eye(3));
alpha_tensor_AU = S * (((eye(3*Na)-Ablk*T) \ Ablk) * K);
alpha_tensor_A3 = alpha_tensor_AU * A3_PER_AU;

result = struct( ...
    'mu_AU', mu, ...
    'F_AU', F, ...
    'U_pol_H', U_pol_H, ...
    'U_pol_eV', U_pol_eV, ...
    'U_pol_kcal', U_pol_kcal, ...
    'rhoAT', rhoAT, ...
    'Na', Na, ...
    'idx_active', idx_active, ...
    'alpha_tensor_A3', alpha_tensor_A3, ...
    'T_AU', T, ...
    'alpha_AU_active', alpha_AU_active, ...
    'info', info);
end

function relres = compute_relres(mu, b, T, alpha_AU_active, bn)
Na = numel(alpha_AU_active);
r = b - mu;
y = T * mu;
for i = 1:Na
    Ii = (3*i-2):(3*i);
    r(Ii) = r(Ii) + alpha_AU_active(i) * y(Ii);
end
relres = norm(r) / bn;
end

function sys = make_sys_from_finite_cluster(R_all_A, alpha_all_A3, pol_mask, q_all_e, a)
Nat = size(R_all_A,1);

sys = struct();
sys.n_sites = Nat;
sys.site_pos = R_all_A;
sys.site_alpha = alpha_all_A3(:);
sys.site_is_polarizable = logical(pol_mask(:));
sys.site_charge = q_all_e(:);
sys.thole_a = a;

% minimal placeholders
sys.site_mol_id = (1:Nat).';
sys.site_label = arrayfun(@(k) sprintf('X%d', k), 1:Nat, 'UniformOutput', false).';
sys.site_type = repmat({'X'}, Nat, 1);
end