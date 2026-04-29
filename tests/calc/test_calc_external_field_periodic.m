function test_calc_external_field_periodic()
%TEST_CALC_EXTERNAL_FIELD_PERIODIC Verify periodic charge external field.
%
% This is an interface, convention, and sanity test, not a hard analytic
% Ewald regression.
%
% It checks:
%   - periodic external field returns n_sites x 3
%   - neutral source-charge requirement is satisfied for +1/-1 pair
%   - target_mask is respected
%   - source_mask is respected
%   - non-neutral periodic source charge errors
%   - periodic k-vectors use the project row-lattice convention
%   - reciprocal/self/surface helper conventions are internally consistent

polsys = local_make_periodic_polsys();

params = struct();
params.use_thole = false;

params.field = struct();
params.field.mode = 'periodic';
params.field.exclude_self = true;
params.field.use_thole_damping = false;
params.field.target_mask = logical(polsys.site_is_polarizable(:));
params.field.source_mask = abs(polsys.site_charge(:)) > 0;
params.field.kspace_mode = 'full';

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

% Direct call should match calc dispatcher and expose diagnostics.
fieldParamsDirect = params.field;
fieldParamsDirect = rmfield(fieldParamsDirect, 'mode');

[E_direct, parts] = thole.induced_field_from_charges_periodic( ...
    polsys, fieldParamsDirect);

assert(norm(E - E_direct, 'fro') < 1e-12, ...
    'calc.compute_external_field should match direct periodic charge-field call.');

assert(parts.nK > 0, ...
    'Periodic external field test should exercise reciprocal k-space.');

assert(isfield(parts, 'kmeta') && isstruct(parts.kmeta), ...
    'Periodic parts should include k-vector metadata.');

assert(isfield(parts.kmeta, 'convention'), ...
    'Periodic k-vector metadata should include convention.');

assert(strcmp(parts.kmeta.convention, 'project_row_H_column_G_HG_2piI'), ...
    'Periodic k-vector enumeration should use the project row-lattice convention.');

assert(isfield(parts, 'lattice') && isstruct(parts.lattice), ...
    'Periodic parts should include canonical lattice metadata.');

assert(norm(parts.lattice.H * parts.lattice.G - 2*pi*eye(3), 'fro') < 1e-12, ...
    'Canonical lattice should satisfy H*G = 2*pi*I.');

assert(norm(parts.lattice.H - polsys.super_lattice, 'fro') < 1e-12, ...
    'Canonical lattice H should match polsys.super_lattice in row convention.');

% Verify returned k-vectors are consistent with k = (G*m).'.
for a = 1:parts.kmeta.num_kvec
    m = parts.kmeta.hkl(a, :).';
    kExpected = (parts.lattice.G * m).';
    assert(norm(parts.kmeta.hkl(a, :)) > 0, ...
        'Returned hkl should exclude k=0.');
    assert(norm(parts.kmeta.knorm(a) - norm(kExpected)) < 1e-12, ...
        'Returned k-vector norm should match G*hkl.');
end

% Check no +/- duplicate hkl representatives.
hkl = parts.kmeta.hkl;
for a = 1:size(hkl, 1)
    hasOpposite = any(all(hkl == -hkl(a, :), 2));
    assert(~hasOpposite, ...
        'Half-space k-vector enumeration should not contain both k and -k.');
end

assert(issorted(parts.kmeta.knorm), ...
    'k-vectors should be sorted by increasing norm.');

assert(all(parts.kmeta.knorm <= params.field.ewald.kcut + 1e-12), ...
    'All k-vectors should lie within kcut.');

% Self/surface helper convention checks folded into this periodic test.
alpha = params.field.ewald.alpha;
Tself = ewald.self_tensor_block_dipole(alpha);
TselfExpected = (4 * alpha^3 / (3 * sqrt(pi))) * eye(3);
assert(norm(Tself - TselfExpected, 'fro') < 1e-14, ...
    'Dipole self block should use the field-operator sign convention.');

TsurfTin = ewald.surface_tensor_block_dipole(parts.lattice, 'tinfoil');
assert(norm(TsurfTin, 'fro') < 1e-14, ...
    'Tinfoil surface block should be zero.');

TsurfVac = ewald.surface_tensor_block_dipole(parts.lattice, 'vacuum');
TsurfVacExpected = -(4 * pi / (3 * parts.lattice.volume)) * eye(3);
assert(norm(TsurfVac - TsurfVacExpected, 'fro') < 1e-14, ...
    'Vacuum surface block should use the field-operator sign convention.');

% Vacuum surface term should be populated and should follow Esurf = -4pi M/3V.
paramsVac = params;
paramsVac.field.ewald.boundary = 'vacuum';
fieldParamsVac = paramsVac.field;
fieldParamsVac = rmfield(fieldParamsVac, 'mode');
[~, partsVac] = thole.induced_field_from_charges_periodic(polsys, fieldParamsVac);

MqExpected = sum(polsys.site_charge(partsVac.source_sites) .* ...
    polsys.site_pos(partsVac.source_sites, :), 1);
EsurfExpected = -(4 * pi / (3 * partsVac.lattice.volume)) * MqExpected;

assert(norm(partsVac.Mq - MqExpected, 2) < 1e-12, ...
    'Vacuum charge surface dipole Mq should match selected source charges.');

assert(norm(partsVac.Esurf_q - EsurfExpected, 2) < 1e-12, ...
    'Vacuum charge surface field should follow Esurf = -4*pi*Mq/(3V).');

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

% Real-only mode should populate only the real contribution and skip k-space.
paramsRealOnly = params;
paramsRealOnly.field.real_only = true;
fieldParamsRealOnly = paramsRealOnly.field;
fieldParamsRealOnly = rmfield(fieldParamsRealOnly, 'mode');

[E_real_only, partsRealOnly] = thole.induced_field_from_charges_periodic( ...
    polsys, fieldParamsRealOnly);

assert(partsRealOnly.real_only, ...
    'real_only diagnostics should record real_only = true.');

assert(partsRealOnly.nK == 0, ...
    'real_only periodic field should skip k-vector enumeration.');

assert(strcmp(partsRealOnly.storage_mode, 'real_only'), ...
    'real_only periodic field should report storage_mode = real_only.');

assert(norm(E_real_only - partsRealOnly.real, 'fro') < 1e-12, ...
    'real_only total field should equal the real-space contribution.');

assert(norm(partsRealOnly.recip, 'fro') < 1e-12, ...
    'real_only reciprocal contribution should be zero.');

assert(norm(partsRealOnly.surf, 'fro') < 1e-12, ...
    'real_only surface contribution should be zero.');

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
    assert(contains(ME.message, 'neutral') || ...
           contains(ME.identifier, 'neutral', 'IgnoreCase', true), ...
        'Non-neutral periodic source error should mention neutrality.');
end

assert(didError, ...
    'Periodic external field should error for non-neutral source charges.');

io.assert_atomic_units(polsys);
end

function polsys = local_make_periodic_polsys()
% Four sites in a mildly triclinic row-lattice cell:
%
%   cart = frac * H
%
% site 1: +1 source
% site 2: -1 source
% site 3: polarizable target
% site 4: neutral spectator
%
% The +1/-1 source mask is neutral, as required by periodic charge Ewald.
%
% The non-orthogonal lattice is intentional: a cubic cell can hide row/column
% convention mistakes because H == H.'.

polsys = struct();

polsys.site_pos = [
     2.0   2.0   2.0
     6.0   2.0   2.0
    10.0   9.0   8.0
    15.0   5.0   4.0
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

% Row-lattice convention: cart = frac * H.
polsys.lattice = [
    20.0   0.0   0.0
     2.0  18.0   0.0
     1.0   3.0  22.0
];

polsys.super_lattice = polsys.lattice;
end