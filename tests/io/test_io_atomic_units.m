function test_io_atomic_units()
%TEST_IO_ATOMIC_UNITS Verify conversion to internal atomic-unit convention.

ANG2BOHR = 1.8897259886;
A3_PER_AU = 0.148184711;

%% Legacy/default behavior: assume Angstrom and Angstrom^3

sys = struct();
sys.site_pos = [
    1.0 0.0 0.0
    0.0 2.0 0.0
];
sys.site_alpha = [
    1.0
    2.0
];
sys.site_charge = [
    1.0
   -1.0
];

sysAU = io.convert_to_atomic_units(sys);

assert(norm(sysAU.site_pos - sys.site_pos * ANG2BOHR, 'fro') < 1e-12, ...
    'Positions were not converted from Angstrom to bohr correctly.');

assert(norm(sysAU.site_alpha - sys.site_alpha / A3_PER_AU) < 1e-12, ...
    'Polarizabilities were not converted from Angstrom^3 to atomic units correctly.');

assert(isfield(sysAU, 'units'), ...
    'Converted system should have a units field.');

assert(strcmp(sysAU.units.length, 'bohr'), ...
    'Length unit should be stamped as bohr.');

assert(strcmp(sysAU.units.alpha, 'atomic_unit'), ...
    'Polarizability unit should be stamped as atomic_unit.');

assert(strcmp(sysAU.units.charge, 'elementary_charge'), ...
    'Charge unit should be stamped as elementary_charge.');

io.assert_atomic_units(sysAU);

%% Already atomic units should be unchanged

sys2 = struct();
sys2.site_pos = [
    3.0 4.0 5.0
    6.0 7.0 8.0
];
sys2.site_alpha = [
    10.0
    20.0
];
sys2.site_charge = [
    0.0
    1.0
];
sys2.units.length = 'bohr';
sys2.units.alpha = 'atomic_unit';
sys2.units.charge = 'elementary_charge';

sys2AU = io.convert_to_atomic_units(sys2);

assert(norm(sys2AU.site_pos - sys2.site_pos, 'fro') < 1e-12, ...
    'Bohr positions should not be changed.');

assert(norm(sys2AU.site_alpha - sys2.site_alpha(:)) < 1e-12, ...
    'Atomic-unit polarizabilities should not be changed.');

io.assert_atomic_units(sys2AU);

%% assert_atomic_units should reject missing metadata

bad = rmfield(sys2AU, 'units');

didFail = false;
try
    io.assert_atomic_units(bad);
catch
    didFail = true;
end

assert(didFail, ...
    'assert_atomic_units should fail when sys.units is missing.');

%% assert_atomic_units should reject wrong metadata

bad = sys2AU;
bad.units.length = 'angstrom';

didFail = false;
try
    io.assert_atomic_units(bad);
catch
    didFail = true;
end

assert(didFail, ...
    'assert_atomic_units should fail for non-bohr length units.');

end