function test_calc_external_field_nonperiodic()
%TEST_CALC_EXTERNAL_FIELD_NONPERIODIC Verify nonperiodic charge field.
%
% We use a tiny atomic-unit system:
%
%   source charge +1 at [0,0,0]
%   target polarizable site at [2,0,0]
%
% With no Thole damping, the electric field at the target is:
%
%   E = q * r / |r|^3 = [2,0,0] / 8 = [0.25,0,0]
%
% where r = r_target - r_source.

polsys = local_make_two_site_polsys();

params = struct();
params.use_thole = false;
params.field = struct();
params.field.mode = 'nonperiodic';
params.field.exclude_self = true;
params.field.use_thole_damping = false;
params.field.target_mask = logical(polsys.site_is_polarizable(:));
params.field.source_mask = abs(polsys.site_charge(:)) > 0;

E = calc.compute_external_field(polsys, params);

assert(isequal(size(E), [2 3]), ...
    'External field should be n_sites x 3.');

expected = [
    0.00 0.00 0.00
    0.25 0.00 0.00
];

assert(norm(E - expected, 'fro') < 1e-12, ...
    'Nonperiodic external field does not match analytic Coulomb field.');

% If target_mask excludes the polarizable target, field should be zero there.
params.field.target_mask = [false; false];

EzeroTarget = calc.compute_external_field(polsys, params);

assert(norm(EzeroTarget, 'fro') < 1e-12, ...
    'External field should be zero when target_mask excludes all sites.');

% If source_mask excludes the charge, field should be zero everywhere.
params.field.target_mask = logical(polsys.site_is_polarizable(:));
params.field.source_mask = [false; false];

EzeroSource = calc.compute_external_field(polsys, params);

assert(norm(EzeroSource, 'fro') < 1e-12, ...
    'External field should be zero when source_mask excludes all sites.');

% Check a negative charge flips the field direction.
polsysNeg = polsys;
polsysNeg.site_charge(1) = -1;

params.field.source_mask = abs(polsysNeg.site_charge(:)) > 0;

Eneg = calc.compute_external_field(polsysNeg, params);

expectedNeg = [
    0.00  0.00 0.00
   -0.25  0.00 0.00
];

assert(norm(Eneg - expectedNeg, 'fro') < 1e-12, ...
    'Negative source charge should flip field direction.');

io.assert_atomic_units(polsys);

end

function polsys = local_make_two_site_polsys()

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
    1.0
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