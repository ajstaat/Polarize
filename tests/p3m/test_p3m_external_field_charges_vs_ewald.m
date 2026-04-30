function test_p3m_external_field_charges_vs_ewald()
%TEST_P3M_EXTERNAL_FIELD_CHARGES_VS_EWALD Compare P3M charge field to Ewald.
%
% This is an early high-level P3M external-field test.
%
% It checks:
%   - p3m.compute_external_field_charges runs on a tiny periodic system
%   - output is finite and zero on non-target rows
%   - mesh charge is conserved
%   - P3M field is broadly aligned with pure Ewald reference
%   - a finer mesh is not worse than a coarser mesh
%
% This is NOT intended to certify production P3M accuracy. It is a smoke /
% regression test while the P3M backend is under construction.

rng(31);

sys = local_make_periodic_charge_system();

targetMask = logical(sys.site_is_polarizable(:));
sourceMask = abs(sys.site_charge(:)) > 0;

assert(abs(sum(sys.site_charge(sourceMask))) < 1e-14, ...
    'Test source charges should be neutral.');
assert(nnz(sourceMask) == 2, ...
    'Test should have two charged source sites.');
assert(nnz(targetMask) >= 3, ...
    'Test should have at least three target sites.');

ewaldParams = struct();
ewaldParams.alpha = 0.30;
ewaldParams.rcut = 9.0;
ewaldParams.kcut = 2.0;
ewaldParams.boundary = 'tinfoil';

% -------------------------------------------------------------------------
% Pure Ewald reference.
% -------------------------------------------------------------------------

fieldEwald = struct();
fieldEwald.ewald = ewaldParams;
fieldEwald.target_mask = targetMask;
fieldEwald.source_mask = sourceMask;
fieldEwald.exclude_self = true;
fieldEwald.use_thole_damping = true;
fieldEwald.real_only = false;
fieldEwald.kspace_mode = 'full';
fieldEwald.k_block_size = 256;
fieldEwald.kspace_memory_limit_gb = 2;
fieldEwald.verbose = false;

[Eewald, partsEwald] = thole.induced_field_from_charges_periodic(sys, fieldEwald);

assert(all(isfinite(Eewald(:))), ...
    'Pure Ewald reference field should be finite.');
assert(norm(Eewald(targetMask, :), 'fro') > 0, ...
    'Pure Ewald reference field should be nonzero on target sites.');
assert(partsEwald.nK > 0, ...
    'Pure Ewald reference should include reciprocal k-vectors.');

% -------------------------------------------------------------------------
% P3M coarse mesh.
% -------------------------------------------------------------------------

optsCoarse = struct();
optsCoarse.ewald = ewaldParams;
optsCoarse.mesh_size = [12 10 8];
optsCoarse.assignment_order = 4;
optsCoarse.target_mask = targetMask;
optsCoarse.source_mask = sourceMask;
optsCoarse.exclude_self = true;
optsCoarse.use_thole_damping = true;
optsCoarse.realspace_backend = 'thole_periodic_real';
optsCoarse.derivative_mode = 'spectral';
optsCoarse.influence_mode = 'ewald';
optsCoarse.deconvolve_assignment = true;
optsCoarse.deconvolution_floor = 1e-8;
optsCoarse.verbose = false;

[Ep3mCoarse, partsCoarse] = p3m.compute_external_field_charges(sys, optsCoarse);

local_assert_p3m_parts_valid(Ep3mCoarse, partsCoarse, sys, targetMask, sourceMask);

% -------------------------------------------------------------------------
% P3M finer mesh.
% -------------------------------------------------------------------------

optsFine = optsCoarse;
optsFine.mesh_size = [18 15 12];

[Ep3mFine, partsFine] = p3m.compute_external_field_charges(sys, optsFine);

local_assert_p3m_parts_valid(Ep3mFine, partsFine, sys, targetMask, sourceMask);

% -------------------------------------------------------------------------
% Compare to Ewald on target rows.
% -------------------------------------------------------------------------

EtEwald = Eewald(targetMask, :);
EtCoarse = Ep3mCoarse(targetMask, :);
EtFine = Ep3mFine(targetMask, :);

relErrCoarse = norm(EtCoarse - EtEwald, 'fro') / max(norm(EtEwald, 'fro'), eps);
relErrFine = norm(EtFine - EtEwald, 'fro') / max(norm(EtEwald, 'fro'), eps);

cosCoarse = local_field_cosine(EtCoarse, EtEwald);
cosFine = local_field_cosine(EtFine, EtEwald);

fprintf('\nP3M external charge-field comparison:\n');
fprintf('  Ewald ||E||_F       = %.12e\n', norm(EtEwald, 'fro'));
fprintf('  coarse ||E||_F      = %.12e | rel err = %.6e | cosine = %.6f\n', ...
    norm(EtCoarse, 'fro'), relErrCoarse, cosCoarse);
fprintf('  fine   ||E||_F      = %.12e | rel err = %.6e | cosine = %.6f\n', ...
    norm(EtFine, 'fro'), relErrFine, cosFine);
fprintf('  coarse mesh         = [%d %d %d]\n', partsCoarse.mesh_size);
fprintf('  fine mesh           = [%d %d %d]\n', partsFine.mesh_size);

% Early P3M smoke tolerances. These are deliberately loose because the
% current influence function is not yet optimized and the mesh is tiny.
assert(cosCoarse > 0.50, ...
    'Coarse P3M field should be broadly aligned with Ewald reference.');
assert(cosFine > 0.65, ...
    'Fine P3M field should be reasonably aligned with Ewald reference.');

assert(relErrFine < 1.25 * relErrCoarse, ...
    'Finer P3M mesh should not be significantly worse than coarse mesh.');

% We also require the fine mesh to be within a broad relative error bound so
% a completely broken normalization cannot slip through.
assert(relErrFine < 1.0, ...
    'Fine P3M external field is too far from Ewald reference: rel err %.3e.', ...
    relErrFine);

% Real-space part should match the pure Ewald real-space part closely
% because compute_external_field_charges reuses the same periodic real cache.
realDiff = norm(partsFine.real - partsEwald.real, 'fro') / ...
    max(norm(partsEwald.real, 'fro'), eps);

assert(realDiff < 1e-12, ...
    'P3M real-space field should match pure Ewald real-space field.');

io.assert_atomic_units(sys);
end

function local_assert_p3m_parts_valid(E, parts, sys, targetMask, sourceMask)
nSites = sys.n_sites;

assert(isequal(size(E), [nSites, 3]), ...
    'P3M field should be N x 3.');

assert(all(isfinite(E(:))), ...
    'P3M field should be finite.');

assert(norm(E(targetMask, :), 'fro') > 0, ...
    'P3M field should be nonzero on target rows.');

assert(norm(E(~targetMask, :), 'fro') == 0, ...
    'P3M field should be exactly zero on non-target rows.');

assert(isfield(parts, 'real') && isequal(size(parts.real), [nSites, 3]), ...
    'parts.real should be N x 3.');
assert(isfield(parts, 'recip') && isequal(size(parts.recip), [nSites, 3]), ...
    'parts.recip should be N x 3.');
assert(isfield(parts, 'surf') && isequal(size(parts.surf), [nSites, 3]), ...
    'parts.surf should be N x 3.');

assert(all(isfinite(parts.real(:))) && ...
       all(isfinite(parts.recip(:))) && ...
       all(isfinite(parts.surf(:))), ...
    'P3M field parts should be finite.');

assert(norm(parts.recip(targetMask, :), 'fro') > 0, ...
    'P3M reciprocal field should be nonzero on target rows.');

assert(abs(parts.qtot - sum(sys.site_charge(sourceMask))) < 1e-14, ...
    'parts.qtot should match selected source charge.');

assert(abs(parts.rho_total - sum(sys.site_charge(sourceMask))) < 1e-13, ...
    'Assigned P3M charge mesh should conserve selected source charge.');

assert(strcmp(parts.lattice_convention, 'project_row_H_column_G_HG_2piI'), ...
    'P3M parts should advertise canonical lattice convention.');

assert(isequal(parts.target_mask(:), targetMask(:)), ...
    'P3M parts target_mask should match input targetMask.');
assert(isequal(parts.source_mask(:), sourceMask(:)), ...
    'P3M parts source_mask should match input sourceMask.');

assert(parts.nK > 0, ...
    'P3M spectral solve should report nonzero reciprocal modes.');
end

function c = local_field_cosine(A, B)
a = A(:);
b = B(:);

den = max(norm(a) * norm(b), eps);
c = dot(a, b) / den;
end

function sys = local_make_periodic_charge_system()
sys = struct();

% Five sites:
%   1-2: fixed charged source sites, nonpolarizable but with physical alpha
%   3-5: polarizable target sites
%
% The cell is mildly triclinic to exercise row-lattice conventions.
sys.site_pos = [
     4.0   3.0   2.0
    13.0   7.0   6.0
     8.0  16.0   7.0
    19.0  10.0  13.0
    23.0  21.0  18.0
];

sys.site_charge = [
    +1.0
    -1.0
     0.0
     0.0
     0.0
];

sys.site_alpha = [
    0.08
    0.07
    0.06
    0.05
    0.04
];

sys.site_is_polarizable = [
    false
    false
    true
    true
    true
];

sys.site_is_active = abs(sys.site_charge) > 0;

sys.n_sites = 5;
sys.thole_a = 0.39;

sys.site_type = {'X'; 'X'; 'X'; 'X'; 'X'};
sys.site_class = {'source'; 'source'; 'pol'; 'pol'; 'pol'};
sys.site_label = {'q+'; 'q-'; 'p1'; 'p2'; 'p3'};
sys.site_mol_id = (1:5).';

sys.units = struct();
sys.units.length = 'bohr';
sys.units.alpha = 'atomic_unit';
sys.units.charge = 'elementary_charge';

sys.is_periodic = true;
sys.periodic_mode = 'periodic';

% Row-lattice convention: cart = frac * H.
sys.lattice = [
    32.0   0.0   0.0
     2.0  31.0   0.0
     1.0   3.0  30.0
];

sys.super_lattice = sys.lattice;
end