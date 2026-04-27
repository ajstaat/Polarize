function test_thole_direct_one_site_energy()
%TEST_THOLE_DIRECT_ONE_SITE_ENERGY Verify direct SCF algebra for one active site.
%
% System:
%   site 1: charged source, nonpolarizable
%   site 2: polarizable target, alpha = 2 au
%
% External field:
%   E(target) = [0.25 0 0]
%
% Dense active-space operator:
%   Tpol = zeros(3,3)
%
% Expected:
%   mu = alpha * E = [0.5 0 0]
%
% Energy:
%   polarization_self      = 0.5 * mu^2 / alpha = +0.0625
%   external_charge_dipole = -mu dot E           = -0.125
%   dipole_dipole          = 0
%   total                  = -0.0625

polsys = local_make_one_site_polsys();

Eext = zeros(polsys.n_sites, 3);
Eext(2, :) = [0.25 0.0 0.0];

scfParams = struct();
scfParams.tol = 1e-12;
scfParams.maxIter = 50;
scfParams.mixing = 1.0;
scfParams.omega = 1.0;
scfParams.verbose = false;

problem = thole.prepare_scf_problem(polsys, Eext, scfParams);

assert(problem.nPolSites == 1, ...
    'Expected exactly one polarizable active site.');

assert(isequal(problem.activeSites(:), 2), ...
    'The only active polarizable site should be site 2.');

assert(isequal(size(problem.Eext_pol_vec), [3 1]), ...
    'Active-space Eext vector should be 3 x 1.');

assert(norm(problem.Eext_pol_vec - [0.25; 0; 0]) < 1e-14, ...
    'Active-space Eext vector is incorrect.');

assert(norm(problem.alpha_pol_vec - 2.0 * ones(3,1)) < 1e-14, ...
    'Active-space alpha vector should repeat alpha for x/y/z.');

% No dipole-dipole coupling in the one-active-site test.
Tpol = zeros(3, 3);

[mu, info] = thole.solve_scf_direct(problem, Tpol);

expectedMu = zeros(polsys.n_sites, 3);
expectedMu(2, :) = [0.5 0.0 0.0];

assert(isequal(size(mu), [polsys.n_sites 3]), ...
    'solve_scf_direct should return full-system N x 3 induced dipoles.');

assert(norm(mu - expectedMu, 'fro') < 1e-12, ...
    'Direct one-site SCF dipole does not match mu = alpha * E.');

assert(isfield(info, 'relres'), ...
    'Direct solver info should contain relres.');

assert(info.relres < 1e-12, ...
    'Direct solver residual should be near zero.');

relres = thole.compute_active_space_relres(problem, Tpol, mu);

assert(relres < 1e-12, ...
    'compute_active_space_relres should report near-zero residual.');

energy = calc.compute_total_energy_active_space(polsys, problem, mu, Eext, Tpol);

assert(abs(energy.polarization_self - 0.0625) < 1e-12, ...
    'polarization_self energy is incorrect.');

assert(abs(energy.external_charge_dipole - (-0.125)) < 1e-12, ...
    'external_charge_dipole energy is incorrect.');

assert(abs(energy.dipole_dipole - 0.0) < 1e-12, ...
    'dipole_dipole energy should be zero when Tpol = 0.');

assert(abs(energy.total - (-0.0625)) < 1e-12, ...
    'Total polarization energy is incorrect.');

assert(isfield(energy, 'total_stationary'), ...
    'Energy output should contain total_stationary.');

assert(isfield(energy, 'stationary_consistency'), ...
    'Energy output should contain stationary_consistency.');

assert(abs(energy.stationary_consistency) < 1e-12, ...
    'Stationary-form energy should match total energy.');

io.assert_atomic_units(polsys);

end

function polsys = local_make_one_site_polsys()

polsys = struct();

polsys.site_pos = [
    0.0 0.0 0.0
    2.0 0.0 0.0
];

polsys.site_charge = [
    +1.0
     0.0
];

polsys.site_alpha = [
    0.0
    2.0
];

polsys.site_is_polarizable = [
    false
    true
];

polsys.site_type = {'X'; 'X'};
polsys.site_class = {'source'; 'target'};
polsys.site_label = {'q1'; 'p1'};
polsys.site_mol_id = [1; 2];
polsys.site_is_active = [true; false];

polsys.n_sites = 2;

polsys.units = struct();
polsys.units.length = 'bohr';
polsys.units.alpha = 'atomic_unit';
polsys.units.charge = 'elementary_charge';

polsys.thole_a = 0.39;

polsys.is_periodic = false;
polsys.periodic_mode = 'nonperiodic';

end