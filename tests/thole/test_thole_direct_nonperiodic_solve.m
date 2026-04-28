function test_thole_direct_nonperiodic_solve()
%TEST_THOLE_NONPERIODIC_DIRECT_WORKFLOW End-to-end finite direct workflow.
%
% This test exercises the restored nonperiodic direct path:
%
%   polsys
%   -> calc.compute_external_field(..., mode='nonperiodic')
%   -> thole.prepare_scf_problem
%   -> thole.assemble_nonperiodic_interaction_matrix
%   -> thole.solve_scf_direct
%   -> calc.compute_total_energy_active_space
%
% System:
%   site 1: +1 charge source, nonpolarizable
%   site 2: -1 charge source, nonpolarizable
%   site 3: polarizable target, alpha = 1.2
%   site 4: polarizable target, alpha = 1.5

polsys = local_make_nonperiodic_workflow_polsys();

fieldParams = struct();
fieldParams.use_thole = false;
fieldParams.field = struct();
fieldParams.field.mode = 'nonperiodic';
fieldParams.field.exclude_self = true;
fieldParams.field.use_thole_damping = false;
fieldParams.field.target_mask = logical(polsys.site_is_polarizable(:));
fieldParams.field.source_mask = abs(polsys.site_charge(:)) > 0;

Eext = calc.compute_external_field(polsys, fieldParams);

assert(isequal(size(Eext), [polsys.n_sites 3]), ...
    'External field should be n_sites x 3.');

assert(all(isfinite(Eext(:))), ...
    'External field should be finite.');

assert(norm(Eext(polsys.site_is_polarizable, :), 'fro') > 0, ...
    'External field on polarizable sites should be nonzero.');

assert(norm(Eext(~polsys.site_is_polarizable, :), 'fro') < 1e-12, ...
    'External field should be zero on non-target source sites.');

scfParams = struct();
scfParams.use_thole = true;
scfParams.softening = 0.0;
scfParams.rcut = Inf;
scfParams.tol = 1e-12;
scfParams.maxIter = 50;
scfParams.mixing = 1.0;
scfParams.omega = 1.0;
scfParams.verbose = false;

problem = thole.prepare_scf_problem(polsys, Eext, scfParams);

assert(problem.nPolSites == 2, ...
    'Expected two polarizable active sites.');

assert(isequal(problem.activeSites(:), [3; 4]), ...
    'Expected polarizable active sites to be sites 3 and 4.');

[Tpol, opinfo] = thole.assemble_nonperiodic_interaction_matrix( ...
    polsys, problem, scfParams);

assert(isequal(size(Tpol), [6 6]), ...
    'Two active polarizable sites should produce a 6 x 6 Tpol.');

assert(norm(Tpol - Tpol.', 'fro') < 1e-12, ...
    'Dense nonperiodic Tpol should be symmetric.');

assert(opinfo.nPolSites == 2, ...
    'Operator info should report two polarizable sites.');

assert(opinfo.nPairBlocks == 1, ...
    'Two polarizable sites should have one pair block.');

assert(opinfo.nPairBlocksKept == 1, ...
    'Pair block should be kept with rcut = Inf.');

assert(opinfo.nPairBlocksSkippedCutoff == 0, ...
    'No pair block should be skipped with rcut = Inf.');

[mu, info] = thole.solve_scf_direct(problem, Tpol);

assert(isequal(size(mu), [polsys.n_sites 3]), ...
    'Direct solver should return full-system n_sites x 3 induced dipoles.');

assert(info.relres < 1e-11, ...
    'Direct solver residual should be small.');

relres = thole.compute_active_space_relres(problem, Tpol, mu);

assert(relres < 1e-11, ...
    'Active-space residual should be small.');

assert(norm(mu(~polsys.site_is_polarizable, :), 'fro') < 1e-12, ...
    'Nonpolarizable charged sites should have zero induced dipoles.');

assert(norm(mu(polsys.site_is_polarizable, :), 'fro') > 0, ...
    'Polarizable sites should have nonzero induced dipoles.');

energy = calc.compute_total_energy_active_space(polsys, problem, mu, Eext, Tpol);

assert(isfinite(energy.total), ...
    'Energy total should be finite.');

assert(abs(energy.stationary_consistency) < 1e-11, ...
    'Stationary energy consistency should be small.');

% Check cutoff behavior separately.
scfParamsCut = scfParams;
scfParamsCut.rcut = 1.0;  % smaller than the site 3--4 separation

[Tcut, opinfoCut] = thole.assemble_nonperiodic_interaction_matrix( ...
    polsys, problem, scfParamsCut);

assert(norm(Tcut, 'fro') < 1e-12, ...
    'Tpol should be zero when rcut excludes the only polarizable pair.');

assert(opinfoCut.use_cutoff, ...
    'Operator info should record that a finite cutoff was used.');

assert(opinfoCut.nPairBlocksKept == 0, ...
    'No pair blocks should be kept with short cutoff.');

assert(opinfoCut.nPairBlocksSkippedCutoff == 1, ...
    'The only pair block should be skipped by cutoff.');

io.assert_atomic_units(polsys);

end

function polsys = local_make_nonperiodic_workflow_polsys()

polsys = struct();

polsys.site_pos = [
    0.0  0.0 0.0   % + charge
    4.0  0.0 0.0   % - charge
    1.5  2.0 0.0   % polarizable
    3.0  2.5 0.0   % polarizable
];

polsys.site_charge = [
    +1.0
    -1.0
     0.0
     0.0
];

polsys.site_alpha = [
    0.0
    0.0
    1.2
    1.5
];

polsys.site_is_polarizable = [
    false
    false
    true
    true
];

polsys.site_type = {'X'; 'X'; 'X'; 'X'};
polsys.site_class = {'qplus'; 'qminus'; 'pol'; 'pol'};
polsys.site_label = {'q+'; 'q-'; 'p1'; 'p2'};
polsys.site_mol_id = [1; 2; 3; 4];
polsys.site_is_active = [true; true; false; false];

polsys.n_sites = 4;

polsys.units = struct();
polsys.units.length = 'bohr';
polsys.units.alpha = 'atomic_unit';
polsys.units.charge = 'elementary_charge';

polsys.thole_a = 0.39;

polsys.is_periodic = false;
polsys.periodic_mode = 'nonperiodic';

end