function test_thole_periodic_ewald_full_vs_blocked_apply()
%TEST_THOLE_PERIODIC_EWALD_FULL_VS_BLOCKED_APPLY Full/blocked k-space consistency.
%
% Checks:
%   1. full and blocked k-space storage produce the same op.apply
%   2. full and blocked operators produce the same SCF solution/energy
%   3. blocked storage is actually selected when requested
%
% This protects the memory-safe blocked path that will matter for larger
% periodic Ewald reference calculations.

rng(12);

sys = local_make_periodic_polsys();
Eext = local_make_external_field(sys);

scfParams = struct();
scfParams.tol = 1e-11;
scfParams.maxIter = 300;
scfParams.mixing = 0.5;
scfParams.omega = 1.0;
scfParams.verbose = false;

problem = thole.prepare_scf_problem(sys, Eext, scfParams);

commonArgs = { ...
    'Mode', 'periodic_ewald', ...
    'Solver', 'gmres', ...
    'Backend', 'auto', ...
    'UseThole', true, ...
    'Rcut', 9.0, ...
    'Alpha', 0.30, ...
    'Kcut', 1.50, ...
    'Boundary', 'tinfoil', ...
    'UseMex', false, ...
    'Profile', false, ...
    'Verbose', false};

opFull = thole.make_polarization_operator( ...
    sys, problem, commonArgs{:}, ...
    'KspaceMode', 'full');

opBlocked = thole.make_polarization_operator( ...
    sys, problem, commonArgs{:}, ...
    'KspaceMode', 'blocked', ...
    'KBlockSize', 2);

assert(strcmp(opFull.mode, 'periodic_ewald'), ...
    'Expected full operator to use periodic_ewald mode.');
assert(strcmp(opBlocked.mode, 'periodic_ewald'), ...
    'Expected blocked operator to use periodic_ewald mode.');

assert(strcmp(opFull.k_cache.storage_mode, 'full'), ...
    'Full operator should use full k-space storage.');
assert(strcmp(opBlocked.k_cache.storage_mode, 'blocked'), ...
    'Blocked operator should use blocked k-space storage.');

assert(opFull.k_cache.num_kvec == opBlocked.k_cache.num_kvec, ...
    'Full and blocked operators should enumerate the same k-vectors.');

assert(opFull.k_cache.num_kvec > 0, ...
    'Test should exercise reciprocal k-space.');

assert(strcmp(opFull.k_cache.convention, 'project_row_H_column_G_HG_2piI'), ...
    'Full operator should use canonical k-vector convention.');
assert(strcmp(opBlocked.k_cache.convention, 'project_row_H_column_G_HG_2piI'), ...
    'Blocked operator should use canonical k-vector convention.');

nDim = 3 * problem.nPolSites;

for trial = 1:5
    muVec = randn(nDim, 1);

    yFull = opFull.apply(muVec);
    yBlocked = opBlocked.apply(muVec);

    relErr = norm(yFull - yBlocked) / max(1, norm(yFull));

    assert(relErr < 1e-11, ...
        'Full and blocked periodic apply disagree on trial %d: rel error %.3e.', ...
        trial, relErr);
end

% Compare SCF solutions and energies.
gmresOpts = struct();
gmresOpts.tol = 1e-11;
gmresOpts.max_iter = 100;
gmresOpts.restart = [];
gmresOpts.verbose = false;

[muFull, infoFull] = thole.solve_scf_gmres(problem, opFull, gmresOpts);
[muBlocked, infoBlocked] = thole.solve_scf_gmres(problem, opBlocked, gmresOpts);

assert(infoFull.converged, ...
    'Full periodic GMRES solve did not converge.');
assert(infoBlocked.converged, ...
    'Blocked periodic GMRES solve did not converge.');

relMuErr = norm(muFull - muBlocked, 'fro') / max(1, norm(muFull, 'fro'));
assert(relMuErr < 1e-10, ...
    'Full and blocked periodic SCF dipoles disagree: rel error %.3e.', relMuErr);

energyFull = calc.compute_total_energy_active_space(sys, problem, muFull, Eext, opFull);
energyBlocked = calc.compute_total_energy_active_space(sys, problem, muBlocked, Eext, opBlocked);

energyDiff = abs(energyFull.total - energyBlocked.total);
assert(energyDiff < 1e-12, ...
    'Full and blocked periodic SCF energies disagree: abs diff %.3e Ha.', energyDiff);

assert(abs(energyFull.stationary_consistency) < 1e-10, ...
    'Full periodic energy stationary consistency is poor.');
assert(abs(energyBlocked.stationary_consistency) < 1e-10, ...
    'Blocked periodic energy stationary consistency is poor.');

io.assert_atomic_units(sys);
end

function sys = local_make_periodic_polsys()
sys = struct();

sys.site_pos = [
     3.0   2.0   2.0
    10.0   6.0   5.0
    17.0  11.0   9.0
     8.0  19.0  14.0
     6.0  15.0  12.0
];

sys.site_charge = zeros(5, 1);

sys.site_alpha = [
    0.08
    0.07
    0.06
    0.05
    0.00
];

sys.site_is_polarizable = [
    true
    true
    true
    true
    false
];

sys.n_sites = 5;
sys.thole_a = 0.39;

sys.site_type = {'X'; 'X'; 'X'; 'X'; 'X'};
sys.site_class = {'pol'; 'pol'; 'pol'; 'pol'; 'spectator'};
sys.site_label = {'p1'; 'p2'; 'p3'; 'p4'; 'n'};
sys.site_mol_id = [1; 2; 3; 4; 5];
sys.site_is_active = sys.site_is_polarizable;

sys.units = struct();
sys.units.length = 'bohr';
sys.units.alpha = 'atomic_unit';
sys.units.charge = 'elementary_charge';

sys.is_periodic = true;
sys.periodic_mode = 'periodic';

% Row-lattice convention: cart = frac * H.
sys.lattice = [
    30.0   0.0   0.0
     2.0  29.0   0.0
     1.0   3.0  28.0
];

sys.super_lattice = sys.lattice;
end

function Eext = local_make_external_field(sys)
Eext = zeros(sys.n_sites, 3);

Eext(1, :) = [ 1.0e-3, -0.5e-3,  0.2e-3];
Eext(2, :) = [-0.4e-3,  0.8e-3, -0.1e-3];
Eext(3, :) = [ 0.3e-3,  0.1e-3,  0.5e-3];
Eext(4, :) = [-0.2e-3, -0.3e-3,  0.4e-3];
end