function test_builder_apply_file_charges()
%TEST_BUILDER_APPLY_FILE_CHARGES Verify oriented frame-based charge mapping.
%
% This test covers Builder-3b:
%   - match_molecule_atoms_by_frame uses an orientable frame convention
%   - apply_molecule_charges(..., 'Mode', 'file', 'Template', ...)
%     maps nonuniform template charges onto a selected complete molecule
%   - template atom order can differ from target atom order
%   - template charges are rescaled to requested total molecular charge
%   - charged molecule is marked active and removed from polarizable mask

sys = local_make_charge_test_sys();

completeIDs = builder.complete_molecule_ids(sys);
assert(numel(completeIDs) >= 2, ...
    'File-charge test system should have at least two complete molecules.');

templateMolID = completeIDs(1);
targetMolID = completeIDs(2);

templateIdxCanonical = builder.site_indices_for_molecule(sys, templateMolID);
targetIdx = builder.site_indices_for_molecule(sys, targetMolID);

% Canonical nonuniform charges in the molecule's site order.
qCanonical = [
    +0.20
    +0.10
    +0.05
    -0.02
    -0.03
    +0.04
    +0.08
    +0.06
    +0.02
];

assert(abs(sum(qCanonical) - 0.50) < 1e-12, ...
    'Template charges in this test should sum to +0.5 before rescaling.');

% Scramble template atom order so the test catches wrong mapping direction.
perm = [3 1 2 6 4 5 9 7 8];

template = struct();
template.site_pos = sys.site_pos(templateIdxCanonical(perm), :);
template.site_type = sys.site_type(templateIdxCanonical(perm));
template.site_charge = qCanonical(perm);

target = struct();
target.site_pos = sys.site_pos(targetIdx, :);
target.site_type = sys.site_type(targetIdx);

distanceTol = 1e-6;  % bohr; translated identical geometries should be exact

referenceAxis = [0 0 1];
primaryAxis = [1 0 0];

map = builder.match_molecule_atoms_by_frame(template, target, ...
    'DistanceTol', distanceTol, ...
    'ReferenceAxis', referenceAxis, ...
    'PrimaryAxis', primaryAxis, ...
    'AmbiguityTol', 1e-12);

assert(numel(map.template_to_target) == numel(templateIdxCanonical), ...
    'Matcher should return one target index per template atom.');

assert(numel(unique(map.template_to_target)) == numel(templateIdxCanonical), ...
    'Matcher target assignment should be one-to-one.');

assert(map.max_distance < distanceTol, ...
    'Frame-matched translated molecules should match nearly exactly.');

expectedTargetCharges = template.site_charge(map.target_to_template);
expectedTargetCharges = expectedTargetCharges * (1.0 / sum(expectedTargetCharges));

% The target molecule has the same canonical site order as qCanonical, so
% after mapping and rescaling we should recover qCanonical scaled to +1.
expectedCanonical = qCanonical * (1.0 / sum(qCanonical));

assert(all(abs(expectedTargetCharges(:) - expectedCanonical(:)) < 1e-12), ...
    'Frame mapping should recover the canonical oriented charge order.');

sys = builder.apply_molecule_charges(sys, targetMolID, ...
    'Mode', 'file', ...
    'Template', template, ...
    'TotalCharges', +1.0, ...
    'SetActive', true, ...
    'DisablePolarizabilityOnCharged', true, ...
    'ZeroExistingCharges', true, ...
    'RequireComplete', true, ...
    'DistanceTol', distanceTol, ...
    'ReferenceAxis', referenceAxis, ...
    'PrimaryAxis', primaryAxis, ...
    'AmbiguityTol', 1e-12, ...
    'Verbose', false);

targetIdx = builder.site_indices_for_molecule(sys, targetMolID);

assert(abs(sum(sys.site_charge(targetIdx)) - 1.0) < 1e-12, ...
    'Target molecule total charge should be rescaled to +1.');

assert(all(abs(sys.site_charge(targetIdx) - expectedCanonical(:)) < 1e-12), ...
    'Target molecule site charges should match oriented frame-mapped template charges.');

unchargedMask = ~ismember(sys.site_mol_id, targetMolID);

assert(all(sys.site_charge(unchargedMask) == 0), ...
    'Non-target molecules should remain uncharged.');

assert(all(sys.site_is_active(targetIdx)), ...
    'Target molecule sites should be marked active.');

assert(~any(sys.site_is_active(unchargedMask)), ...
    'Non-target molecule sites should not be marked active.');

assert(~any(sys.site_is_polarizable(targetIdx)), ...
    'Target molecule sites should be removed from polarizable mask.');

assert(all(sys.site_alpha(targetIdx) == 0), ...
    'Target molecule polarizabilities should be zeroed.');

assert(all(sys.site_is_polarizable(unchargedMask)), ...
    'Non-target molecule sites should remain polarizable.');

assert(all(sys.site_alpha(unchargedMask) > 0), ...
    'Non-target molecule polarizabilities should remain positive.');

assert(isequal(sys.charged_molecules(:), targetMolID), ...
    'charged_molecules bookkeeping should record target molecule.');

assert(isequal(sys.charged_molecule_total_charges(:), 1.0), ...
    'charged_molecule_total_charges should record requested total charge.');

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

filename = fullfile(tmpDir, 'POSCAR_two_propene_file_charge_test');

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

fprintf(fid, 'Two propene-like molecules for file charge test\n');
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