function test_thole_induced_field_from_dipoles_nonperiodic()
%TEST_THOLE_INDUCED_FIELD_FROM_DIPOLES_NONPERIODIC Verify nonperiodic dipole field.
%
% This anchors the nonperiodic dipole-field convention used by:
%
%   thole.induced_field_from_dipoles_thole
%   thole.dipole_tensor_block
%
% Important:
%   induced_field_from_dipoles_thole currently always applies Thole damping.
%   Therefore the main field test compares against dipole_tensor_block with
%   use_thole = true, rather than against the bare point-dipole analytic
%   value.
%
% We also test the bare no-Thole tensor separately through dipole_tensor_block.

sys = local_make_two_site_sys();

mu = zeros(sys.n_sites, 3);
mu(1, :) = [1.0 0.0 0.0];

params = struct();
params.softening = 0.0;
params.target_mask = [false; true];
params.source_mask = [true; false];

E = thole.induced_field_from_dipoles_thole(sys, mu, params);

assert(isequal(size(E), [sys.n_sites 3]), ...
    'Induced field should be n_sites x 3.');

% Expected damped field from the same pair tensor convention.
opts = struct();
opts.softening = 0.0;
opts.use_thole = true;

T21 = thole.dipole_tensor_block( ...
    sys.site_pos(2, :), ...
    sys.site_pos(1, :), ...
    sys.site_alpha(2), ...
    sys.site_alpha(1), ...
    sys.thole_a, ...
    opts);

expected = zeros(sys.n_sites, 3);
expected(2, :) = (T21 * [1; 0; 0]).';

assert(norm(E - expected, 'fro') < 1e-12, ...
    'Dipole field should match dipole_tensor_block for the same Thole-damped pair.');

% Flip dipole direction.
mu(1, :) = [-1.0 0.0 0.0];

Eflip = thole.induced_field_from_dipoles_thole(sys, mu, params);

assert(norm(Eflip + expected, 'fro') < 1e-12, ...
    'Dipole-field response should be linear in mu.');

% Excluding source should zero the field.
paramsNoSource = params;
paramsNoSource.source_mask = [false; false];

EzeroSource = thole.induced_field_from_dipoles_thole(sys, mu, paramsNoSource);

assert(norm(EzeroSource, 'fro') < 1e-12, ...
    'Field should be zero when source_mask excludes all dipoles.');

% Excluding target should zero the field.
paramsNoTarget = params;
paramsNoTarget.target_mask = [false; false];

EzeroTarget = thole.induced_field_from_dipoles_thole(sys, mu, paramsNoTarget);

assert(norm(EzeroTarget, 'fro') < 1e-12, ...
    'Field should be zero when target_mask excludes all targets.');

% Direct bare point-dipole tensor check.
%
% Geometry:
%   source at [0,0,0]
%   target at [2,0,0]
%   mu = [1,0,0]
%
% Bare field:
%   E = (3 rhat rhat' - I) mu / r^3
%     = [0.25, 0, 0]
optsNoThole = struct();
optsNoThole.softening = 0.0;
optsNoThole.use_thole = false;

Tbare = thole.dipole_tensor_block( ...
    sys.site_pos(2, :), ...
    sys.site_pos(1, :), ...
    sys.site_alpha(2), ...
    sys.site_alpha(1), ...
    sys.thole_a, ...
    optsNoThole);

assert(norm(Tbare * [1; 0; 0] - [0.25; 0; 0]) < 1e-12, ...
    'dipole_tensor_block no-damping tensor is inconsistent with analytic field.');

io.assert_atomic_units(sys);

end

function sys = local_make_two_site_sys()

sys = struct();

sys.site_pos = [
    0.0 0.0 0.0
    2.0 0.0 0.0
];

sys.site_charge = [0.0; 0.0];

sys.site_alpha = [1.0; 1.0];

sys.site_is_polarizable = [true; true];

sys.site_type = {'X'; 'X'};
sys.site_class = {'source'; 'target'};
sys.site_label = {'d1'; 'p1'};
sys.site_mol_id = [1; 2];
sys.site_is_active = [true; true];

sys.n_sites = 2;

sys.units = struct();
sys.units.length = 'bohr';
sys.units.alpha = 'atomic_unit';
sys.units.charge = 'elementary_charge';

sys.thole_a = 0.39;

sys.is_periodic = false;
sys.periodic_mode = 'nonperiodic';

end