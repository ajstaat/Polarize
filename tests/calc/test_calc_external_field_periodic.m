function test_calc_external_field_periodic()
%TEST_CALC_EXTERNAL_FIELD_PERIODIC Verify periodic charge external field.
%
% This is an interface/sanity test, not a hard analytic Ewald regression.
%
% It checks:
%   - periodic external field returns n_sites x 3
%   - neutral source-charge requirement is satisfied for +1/-1 pair
%   - target_mask is respected
%   - source_mask is respected
%   - non-neutral periodic source charge errors

polsys = local_make_periodic_polsys();

params = struct();
params.use_thole = false;
params.field = struct();
params.field.mode = 'periodic';
params.field.exclude_self = true;
params.field.use_thole_damping = false;
params.field.target_mask = logical(polsys.site_is_polarizable(:));
params.field.source_mask = abs(polsys.site_charge(:)) > 0;

params.field.ewald = struct();
params.field.ewald.alpha = 0.35;
params.field.ewald.rcut = 8.0;
params.field.ewald.kcut = 3.0;
params.field.ewald.boundary = 'tinfoil';

E = calc.compute_external_field(polsys, params);

assert(isequal(size(E), [polsys.n_sites 3]), ...
    'Periodic external field should be n_sites x 3.');

assert(all(isfinite(E(:))), ...
    'Periodic external field should be finite.');

targetIdx = find(polsys.site_is_polarizable);

assert(norm(E(targetIdx, :), 'fro') > 0, ...
    'Periodic external field at the polarizable target should be nonzero.');

sourceIdx = find(abs(polsys.site_charge) > 0);

assert(all(vecnorm(E(sourceIdx, :), 2, 2) == 0), ...
    'Periodic external field should be zero on non-target source sites.');

% target_mask excludes all sites -> zero field everywhere.
paramsNoTarget = params;
paramsNoTarget.field.target_mask = false(polsys.n_sites, 1);

E_no_target = calc.compute_external_field(polsys, paramsNoTarget);

assert(norm(E_no_target, 'fro') < 1e-12, ...
    'Periodic external field should be zero when target_mask excludes all sites.');

% source_mask excludes charges -> zero field everywhere.
paramsNoSource = params;
paramsNoSource.field.source_mask = false(polsys.n_sites, 1);

E_no_source = calc.compute_external_field(polsys, paramsNoSource);

assert(norm(E_no_source, 'fro') < 1e-12, ...
    'Periodic external field should be zero when source_mask excludes all sites.');

% Non-neutral periodic source charges should error.
polsysBad = polsys;
polsysBad.site_charge(2) = 0.0;  % source mask now has net +1 if recomputed

paramsBad = params;
paramsBad.field.source_mask = abs(polsysBad.site_charge(:)) > 0;

didError = false;

try
    calc.compute_external_field(polsysBad, paramsBad);
catch ME
    didError = true;
    assert(contains(ME.message, 'neutral') || contains(ME.identifier, 'neutral', 'IgnoreCase', true), ...
        'Non-neutral periodic source error should mention neutrality.');
end

assert(didError, ...
    'Periodic external field should error for non-neutral source charges.');

io.assert_atomic_units(polsys);

end

function polsys = local_make_periodic_polsys()

% Four sites in a simple cubic 20 bohr cell:
%   site 1: +1 source
%   site 2: -1 source
%   site 3: polarizable target
%   site 4: neutral spectator
%
% The +1/-1 source mask is neutral, as required by periodic charge Ewald.

polsys = struct();

polsys.site_pos = [
     2.0  2.0  2.0
     6.0  2.0  2.0
    10.0  9.0  8.0
    15.0  5.0  4.0
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
    1.0
    1.0
];

polsys.site_is_polarizable = [
    false
    false
    true
    false
];

polsys.site_type = {'X'; 'X'; 'X'; 'X'};
polsys.site_class = {'source'; 'source'; 'target'; 'spectator'};
polsys.site_label = {'q+'; 'q-'; 'p'; 'n'};

polsys.site_mol_id = [1; 2; 3; 4];
polsys.site_is_active = [true; true; false; false];

polsys.n_sites = 4;

polsys.units = struct();
polsys.units.length = 'bohr';
polsys.units.alpha = 'atomic_unit';
polsys.units.charge = 'elementary_charge';

polsys.thole_a = 0.39;

polsys.is_periodic = true;
polsys.periodic_mode = 'periodic';

polsys.lattice = 20.0 * eye(3);
polsys.super_lattice = polsys.lattice;

end