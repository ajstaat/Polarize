function test_thole_direct_two_site_coupled()
%TEST_THOLE_DIRECT_TWO_SITE_COUPLED Verify dense direct SCF with dipole coupling.
%
% This test connects the pieces restored so far:
%
%   prepare_scf_problem
%   dipole_tensor_block
%   solve_scf_direct
%   compute_active_space_relres
%   compute_total_energy_active_space
%
% System:
%   site 1: polarizable, alpha = 1.0
%   site 2: polarizable, alpha = 1.5
%
% External field:
%   E1 = [0.10 0 0]
%   E2 = [0.05 0 0]
%
% Coupling:
%   T12/T21 from thole.dipole_tensor_block
%
% Expected direct solve:
%   (I - A T) mu = A E
%
% The test checks:
%   - direct solver agrees with an independently built active-space solve
%   - residual is small
%   - energy stationary consistency is small
%   - induced dipoles are nonzero on both active sites

polsys = local_make_two_site_polsys();

Eext = [
    0.10 0.0 0.0
    0.05 0.0 0.0
];

scfParams = struct();
scfParams.tol = 1e-12;
scfParams.maxIter = 50;
scfParams.mixing = 1.0;
scfParams.omega = 1.0;
scfParams.verbose = false;

problem = thole.prepare_scf_problem(polsys, Eext, scfParams);

assert(problem.nPolSites == 2, ...
    'Expected two polarizable active sites.');

assert(isequal(problem.activeSites(:), [1; 2]), ...
    'Expected both sites to be active polarizable sites.');

Tpol = local_build_dense_Tpol(polsys);

[mu, info] = thole.solve_scf_direct(problem, Tpol);

assert(isequal(size(mu), [polsys.n_sites 3]), ...
    'solve_scf_direct should return full-system N x 3 induced dipoles.');

assert(isfield(info, 'relres'), ...
    'Direct solver info should contain relres.');

assert(info.relres < 1e-11, ...
    'Direct solver residual should be small.');

relres = thole.compute_active_space_relres(problem, Tpol, mu);

assert(relres < 1e-11, ...
    'compute_active_space_relres should report small residual.');

muExpected = local_independent_direct_solution(problem, Tpol);
muExpectedFull = zeros(polsys.n_sites, 3);
muExpectedFull(problem.activeSites, :) = util.unstack_xyz(muExpected);

assert(norm(mu - muExpectedFull, 'fro') < 1e-12, ...
    'Direct SCF solution does not match independent active-space solve.');

assert(norm(mu(1, :)) > 0, ...
    'Site 1 induced dipole should be nonzero.');

assert(norm(mu(2, :)) > 0, ...
    'Site 2 induced dipole should be nonzero.');

energy = calc.compute_total_energy_active_space(polsys, problem, mu, Eext, Tpol);

assert(isfinite(energy.total), ...
    'Energy total should be finite.');

assert(isfield(energy, 'stationary_consistency'), ...
    'Energy output should contain stationary_consistency.');

assert(abs(energy.stationary_consistency) < 1e-11, ...
    'Stationary energy consistency should be small.');

io.assert_atomic_units(polsys);

end

function Tpol = local_build_dense_Tpol(sys)

nPol = nnz(sys.site_is_polarizable);
assert(nPol == 2, 'This helper expects exactly two polarizable sites.');

active = find(sys.site_is_polarizable);
Tpol = zeros(3*nPol, 3*nPol);

opts = struct();
opts.softening = 0.0;
opts.use_thole = true;

for a = 1:nPol
    i = active(a);

    for b = 1:nPol
        j = active(b);

        if i == j
            continue;
        end

        Tij = thole.dipole_tensor_block( ...
            sys.site_pos(i, :), ...
            sys.site_pos(j, :), ...
            sys.site_alpha(i), ...
            sys.site_alpha(j), ...
            sys.thole_a, ...
            opts);

        rows = (3*(a-1)+1):(3*a);
        cols = (3*(b-1)+1):(3*b);

        Tpol(rows, cols) = Tij;
    end
end

end

function muVec = local_independent_direct_solution(problem, Tpol)

A = diag(problem.alpha_pol_vec);

M = eye(numel(problem.alpha_pol_vec)) - A * Tpol;
rhs = A * problem.Eext_pol_vec;

muVec = M \ rhs;

end

function polsys = local_make_two_site_polsys()

polsys = struct();

polsys.site_pos = [
    0.0 0.0 0.0
    3.0 0.0 0.0
];

polsys.site_charge = [0.0; 0.0];

polsys.site_alpha = [1.0; 1.5];

polsys.site_is_polarizable = [true; true];

polsys.site_type = {'X'; 'X'};
polsys.site_class = {'p1'; 'p2'};
polsys.site_label = {'p1'; 'p2'};
polsys.site_mol_id = [1; 2];
polsys.site_is_active = [true; true];

polsys.n_sites = 2;

polsys.units = struct();
polsys.units.length = 'bohr';
polsys.units.alpha = 'atomic_unit';
polsys.units.charge = 'elementary_charge';

polsys.thole_a = 0.39;

polsys.is_periodic = false;
polsys.periodic_mode = 'nonperiodic';

end