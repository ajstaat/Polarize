function test_builder_apply_uniform_charges()
%TEST_BUILDER_APPLY_UNIFORM_CHARGES Verify uniform charge assignment workflow.
%
% This tests the Builder-3a workflow:
%   complete selected pair
%   -> assign +1/-1 uniform charges
%   -> mark selected molecules active
%   -> remove selected molecule sites from induced-dipole active mask
%   -> preserve physical site_alpha values for Thole damping/smearing

sys = local_make_charge_test_sys();
sysBefore = sys;

completeIDs = builder.complete_molecule_ids(sys);
assert(numel(completeIDs) >= 2, ...
    'Charge test system should have at least two complete molecules.');

molIDs = completeIDs(1:2).';
qMol = [+1 -1];

idx1 = builder.site_indices_for_molecule(sys, molIDs(1));
idx2 = builder.site_indices_for_molecule(sys, molIDs(2));

assert(all(sys.site_is_polarizable(idx1)), ...
    'Molecule 1 should start polarizable.');
assert(all(sys.site_is_polarizable(idx2)), ...
    'Molecule 2 should start polarizable.');
assert(all(sys.site_alpha(idx1) > 0), ...
    'Molecule 1 should start with positive physical polarizabilities.');
assert(all(sys.site_alpha(idx2) > 0), ...
    'Molecule 2 should start with positive physical polarizabilities.');
assert(all(sys.site_charge == 0), ...
    'System should start neutral site-by-site.');

sys = builder.apply_molecule_charges(sys, molIDs, ...
    'Mode', 'uniform', ...
    'TotalCharges', qMol, ...
    'SetActive', true, ...
    'DisablePolarizabilityOnCharged', true, ...
    'ZeroExistingCharges', true, ...
    'RequireComplete', true, ...
    'Verbose', false);

idx1 = builder.site_indices_for_molecule(sys, molIDs(1));
idx2 = builder.site_indices_for_molecule(sys, molIDs(2));

assert(abs(sum(sys.site_charge(idx1)) - qMol(1)) < 1e-12, ...
    'First charged molecule should sum to +1.');
assert(abs(sum(sys.site_charge(idx2)) - qMol(2)) < 1e-12, ...
    'Second charged molecule should sum to -1.');

assert(all(abs(sys.site_charge(idx1) - qMol(1) / numel(idx1)) < 1e-12), ...
    'First molecule charge should be distributed uniformly.');
assert(all(abs(sys.site_charge(idx2) - qMol(2) / numel(idx2)) < 1e-12), ...
    'Second molecule charge should be distributed uniformly.');

chargedMask = ismember(sys.site_mol_id, molIDs);

assert(all(sys.site_is_active(chargedMask)), ...
    'Charged molecule sites should be marked active.');
assert(~any(sys.site_is_active(~chargedMask)), ...
    'Uncharged molecule sites should not be marked active.');

assert(isequal(sys.active_molecules(:), molIDs(:)), ...
    'active_molecules should match charged molecule IDs.');
assert(numel(sys.active_site_indices) == 2, ...
    'Expected one active_site_indices entry per charged molecule.');
assert(isequal(sys.active_site_indices{1}(:), idx1(:)), ...
    'First active site-index list should match first charged molecule.');
assert(isequal(sys.active_site_indices{2}(:), idx2(:)), ...
    'Second active site-index list should match second charged molecule.');

assert(~any(sys.site_is_polarizable(chargedMask)), ...
    'Charged molecule sites should be removed from induced-dipole polarizable mask.');

assert(all(sys.site_alpha(chargedMask) > 0), ...
    'Charged molecule physical site_alpha values should be preserved for Thole damping.');

assert(norm(sys.site_alpha(chargedMask) - sysBefore.site_alpha(chargedMask), inf) < 1e-14, ...
    'Charged molecule site_alpha values should be unchanged by charge assignment.');

assert(all(sys.site_is_polarizable(~chargedMask)), ...
    'Uncharged molecule sites should remain polarizable.');
assert(all(sys.site_alpha(~chargedMask) > 0), ...
    'Uncharged molecule site polarizabilities should remain positive.');

assert(abs(sum(sys.site_charge) - sum(qMol)) < 1e-12, ...
    'Total system charge should equal sum of assigned molecular charges.');

assert(isequal(sys.charged_molecules(:), molIDs(:)), ...
    'charged_molecules bookkeeping should match assigned molecule IDs.');
assert(isequal(sys.charged_molecule_total_charges(:), qMol(:)), ...
    'charged_molecule_total_charges bookkeeping should match input charges.');

io.assert_atomic_units(sys);
end

function sys = local_make_charge_test_sys()
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
end

function model = local_make_model()
model = struct();
model.thole_a = 0.39;
model.alpha_units = 'angstrom^3';

model.polarizable_classes = { ...
    'C_deg3' ...
    'C_deg4' ...
    'H_on_C_deg3' ...
    'H_on_C_deg4'};

model.alpha_by_class = struct();
model.alpha_by_class.C_deg3 = 1.750;
model.alpha_by_class.C_deg4 = 1.750;
model.alpha_by_class.H_on_C_deg3 = 0.696;
model.alpha_by_class.H_on_C_deg4 = 0.696;
end

function [filename, tmpDir] = local_write_two_propene_poscar()
tmpDir = tempname;
mkdir(tmpDir);

filename = fullfile(tmpDir, 'POSCAR_two_propene_charge_test');

xyz0 = [
    5.0000 5.0000 5.0000
    6.3400 5.0000 5.0000
    7.0900 6.2990 5.0000
    4.4550 5.9439 5.0000
    4.4550 4.0561 5.0000
    6.8850 4.0561 5.0000
    6.3817 7.1275 5.0000
    7.7166 6.3568 5.8900
    7.7166 6.3568 4.1100
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

fprintf(fid, 'Two propene-like molecules for charge test\n');
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