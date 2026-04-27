function test_builder_extract_polarization_system()
%TEST_BUILDER_EXTRACT_POLARIZATION_SYSTEM Verify solver-facing system extraction.
%
% This tests that charged/active/polarizable fields survive extraction and
% that periodic/nonperiodic mode flags are set consistently.

sys = local_make_charged_sys();

nSites = sys.n_sites;
chargedMask = ismember(sys.site_mol_id, sys.charged_molecules);

polsys = builder.extract_polarization_system(sys, struct( ...
    'mode', 'periodic'));

io.assert_atomic_units(polsys);

assert(polsys.n_sites == nSites, ...
    'Extracted polsys should retain all sites by default.');

assert(isfield(polsys, 'is_periodic') && polsys.is_periodic, ...
    'Periodic extraction should set polsys.is_periodic = true.');

assert(strcmp(polsys.periodic_mode, 'periodic'), ...
    'Periodic extraction should set periodic_mode = periodic.');

assert(isfield(polsys, 'site_pos') && isequal(size(polsys.site_pos), size(sys.site_pos)), ...
    'polsys.site_pos should match sys.site_pos size.');

assert(isequal(polsys.site_charge(:), sys.site_charge(:)), ...
    'polsys.site_charge should preserve system charges.');

assert(isequal(logical(polsys.site_is_polarizable(:)), logical(sys.site_is_polarizable(:))), ...
    'polsys.site_is_polarizable should preserve polarizable mask.');

assert(isequal(logical(polsys.site_is_active(:)), logical(sys.site_is_active(:))), ...
    'polsys.site_is_active should preserve active mask.');

assert(all(abs(polsys.site_charge(chargedMask)) > 0), ...
    'Charged molecule sites should carry nonzero charge in polsys.');

assert(~any(polsys.site_is_polarizable(chargedMask)), ...
    'Charged molecule sites should remain nonpolarizable in polsys.');

assert(all(polsys.site_alpha(chargedMask) == 0), ...
    'Charged molecule site alphas should remain zero in polsys.');

assert(all(polsys.site_alpha(~chargedMask) > 0), ...
    'Uncharged molecule site alphas should remain positive in polsys.');

assert(isfield(polsys, 'lattice') && isequal(polsys.lattice, sys.super_lattice), ...
    'Periodic polsys.lattice should match sys.super_lattice.');

assert(isfield(polsys, 'super_lattice') && isequal(polsys.super_lattice, sys.super_lattice), ...
    'polsys.super_lattice should match sys.super_lattice.');

% Compatibility mode through params.ewald.mode.
polsys2 = builder.extract_polarization_system(sys, struct( ...
    'ewald', struct('mode', 'periodic')));

assert(polsys2.is_periodic, ...
    'params.ewald.mode = periodic should produce periodic polsys.');

% Nonperiodic mode.
polsysNP = builder.extract_polarization_system(sys, struct( ...
    'mode', 'nonperiodic'));

assert(~polsysNP.is_periodic, ...
    'Nonperiodic extraction should set polsys.is_periodic = false.');

assert(strcmp(polsysNP.periodic_mode, 'nonperiodic'), ...
    'Nonperiodic extraction should set periodic_mode = nonperiodic.');

% Active-only extraction should retain only charged active sites.
polsysActive = builder.extract_polarization_system(sys, struct( ...
    'mode', 'periodic', ...
    'active_only', true));

activeMask = logical(sys.site_is_active(:));

assert(polsysActive.n_sites == nnz(activeMask), ...
    'active_only extraction should retain only active sites.');

assert(isequal(polsysActive.site_charge(:), sys.site_charge(activeMask)), ...
    'active_only extraction should preserve active-site charges.');

assert(all(abs(polsysActive.site_charge) > 0), ...
    'All active-only extracted sites should be charged in this test.');

io.assert_atomic_units(polsysActive);

end

function sys = local_make_charged_sys()

[filename, tmpDir] = local_write_two_propene_poscar();
cleanup = onCleanup(@() local_cleanup(tmpDir)); %#ok<NASGU>

crystal = io.import_contcar_as_crystal(filename, ...
    'BondScale', 1.20, ...
    'SortMolecules', false);

model = local_make_model();

sys = builder.make_crystal_system(crystal, model, struct( ...
    'supercell_size', [2 1 1], ...
    'bondScale', 1.20, ...
    'verbose', false));

completeIDs = builder.complete_molecule_ids(sys);
molIDs = completeIDs(1:2).';

sys = builder.apply_molecule_charges(sys, molIDs, ...
    'Mode', 'uniform', ...
    'TotalCharges', [+1 -1], ...
    'SetActive', true, ...
    'DisablePolarizabilityOnCharged', true, ...
    'ZeroExistingCharges', true, ...
    'RequireComplete', true, ...
    'Verbose', false);

end

function model = local_make_model()

model = struct();
model.thole_a = 0.39;
model.alpha_units = 'angstrom^3';

model.polarizable_classes = {
    'C_deg3'
    'C_deg4'
    'H_on_C_deg3'
    'H_on_C_deg4'
};

model.alpha_by_class = struct();
model.alpha_by_class.C_deg3 = 1.750;
model.alpha_by_class.C_deg4 = 1.750;
model.alpha_by_class.H_on_C_deg3 = 0.696;
model.alpha_by_class.H_on_C_deg4 = 0.696;

end

function [filename, tmpDir] = local_write_two_propene_poscar()

tmpDir = tempname;
mkdir(tmpDir);

filename = fullfile(tmpDir, 'POSCAR_two_propene_extract_test');

xyz0 = [
    5.0000  5.0000  5.0000
    6.3400  5.0000  5.0000
    7.0900  6.2990  5.0000

    4.4550  5.9439  5.0000
    4.4550  4.0561  5.0000
    6.8850  4.0561  5.0000

    6.3817  7.1275  5.0000
    7.7166  6.3568  5.8900
    7.7166  6.3568  4.1100
];

xyzA = xyz0;
xyzB = xyz0 + [0 8 0];

allMols = {xyzA, xyzB};

xyzC_all = [];
xyzH_all = [];

for m = 1:numel(allMols)
    xyzM = allMols{m};
    xyzC_all = [xyzC_all; xyzM(1:3, :)]; %#ok<AGROW>
    xyzH_all = [xyzH_all; xyzM(4:9, :)]; %#ok<AGROW>
end

xyz = [xyzC_all; xyzH_all];

fid = fopen(filename, 'w');
assert(fid > 0, 'Failed to open temporary POSCAR.');

fprintf(fid, 'Two propene-like molecules for extract test\n');
fprintf(fid, '1.0\n');
fprintf(fid, '20.0 0.0 0.0\n');
fprintf(fid, '0.0 20.0 0.0\n');
fprintf(fid, '0.0 0.0 12.0\n');
fprintf(fid, 'C H\n');
fprintf(fid, '6 12\n');
fprintf(fid, 'Cartesian\n');

for i = 1:size(xyz, 1)
    fprintf(fid, '%.10f %.10f %.10f\n', xyz(i,1), xyz(i,2), xyz(i,3));
end

fclose(fid);

end

function local_cleanup(tmpDir)

if exist(tmpDir, 'dir')
    rmdir(tmpDir, 's');
end

end