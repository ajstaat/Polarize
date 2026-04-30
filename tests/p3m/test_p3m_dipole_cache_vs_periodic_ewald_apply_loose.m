function test_p3m_dipole_cache_vs_periodic_ewald_apply_loose()
%TEST_P3M_DIPOLE_CACHE_VS_PERIODIC_EWALD_APPLY_LOOSE Compare P3M recip to Ewald recip.
%
% This is an early P3M reciprocal dipole-field validation test.
%
% It compares:
%
%   P3M reciprocal mesh dipole field
%
% against:
%
%   periodic Ewald full operator
% - periodic Ewald no-k operator
%
% so the reference is approximately the pure Ewald reciprocal piece:
%
%   Trecip_Ewald = Tfull_Ewald - Treal_self_surface_Ewald
%
% The comparison is intentionally loose because:
%   - the mesh is tiny
%   - the influence function is not yet optimized
%   - this test is intended to catch sign/normalization/convention failures,
%     not certify production P3M accuracy.

rng(42);

sys = local_make_periodic_polsys();

Eext = zeros(sys.n_sites, 3);

scfParams = struct();
scfParams.tol = 1e-11;
scfParams.maxIter = 200;
scfParams.mixing = 0.5;
scfParams.omega = 1.0;
scfParams.verbose = false;

problem = thole.prepare_scf_problem(sys, Eext, scfParams);

commonArgs = { ...
    'Mode', 'periodic_ewald', ...
    'Solver', 'gmres', ...
    'Backend', 'auto', ...
    'UseThole', false, ...
    'Rcut', 9.0, ...
    'Alpha', 0.30, ...
    'Boundary', 'tinfoil', ...
    'KspaceMode', 'full', ...
    'UseMex', false, ...
    'Profile', false, ...
    'Verbose', false};

opFull = thole.make_polarization_operator(sys, problem, ...
    commonArgs{:}, ...
    'Kcut', 2.0);

% Build an operator with effectively no reciprocal vectors.
%
% Kcut must remain positive for argument validation, so use a tiny positive
% cutoff. The test asserts that no reciprocal k-vectors are actually present.
opNoK = thole.make_polarization_operator(sys, problem, ...
    commonArgs{:}, ...
    'Kcut', 1.0e-12);

assert(strcmp(opFull.mode, 'periodic_ewald'), ...
    'Expected periodic_ewald full reference operator.');
assert(strcmp(opNoK.mode, 'periodic_ewald'), ...
    'Expected periodic_ewald no-k reference operator.');

assert(isfield(opFull.info, 'nK') && opFull.info.nK > 0, ...
    'Full reference operator should have reciprocal k-vectors.');
assert(isfield(opNoK.info, 'nK') && opNoK.info.nK == 0, ...
    'No-k reference operator should have zero reciprocal k-vectors.');

opts = struct();
opts.ewald = struct();
opts.ewald.alpha = 0.30;
opts.mesh_size = [24 22 20];
opts.assignment_order = 4;
opts.target_mask = logical(sys.site_is_polarizable(:));
opts.source_mask = logical(sys.site_is_polarizable(:));
opts.derivative_mode = 'spectral';
opts.influence_mode = 'ewald';
opts.deconvolve_assignment = true;
opts.deconvolution_floor = 1e-8;
opts.verbose = false;

cache = p3m.build_dipole_cache(sys, opts);

nPol = problem.nPolSites;
nSites = sys.n_sites;

% Active-space test dipoles.
muPol = [
     0.10  -0.20   0.30
    -0.40   0.50  -0.60
     0.70   0.10  -0.20
    -0.30  -0.10   0.40
];

assert(size(muPol, 1) == nPol, ...
    'Test dipole matrix should match nPol.');

muVec = util.stack_xyz(muPol);

% Ewald reciprocal-only reference by subtraction.
EfullVec = opFull.apply(muVec);
EnoKVec = opNoK.apply(muVec);
ErecipVec = EfullVec - EnoKVec;

ErecipPol = util.unstack_xyz(ErecipVec);

ErecipFull = zeros(nSites, 3);
ErecipFull(problem.activeSites, :) = ErecipPol;

% P3M reciprocal cache apply.
muFull = zeros(nSites, 3);
muFull(problem.activeSites, :) = muPol;

[Ep3m, partsP3M] = p3m.apply_dipole_cache(cache, muFull);

EtEwald = ErecipFull(problem.activeSites, :);
EtP3M = Ep3m(problem.activeSites, :);

assert(all(isfinite(EtP3M(:))), ...
    'P3M cached reciprocal field should be finite.');
assert(norm(EtP3M, 'fro') > 0, ...
    'P3M cached reciprocal field should be nonzero.');

assert(norm(partsP3M.P_total - sum(muPol, 1)) < 1e-13, ...
    'P3M cached reciprocal scatter should conserve total dipole.');

assert(norm(EtEwald, 'fro') > 0, ...
    'Ewald reciprocal-by-subtraction reference should be nonzero.');

relErr = norm(EtP3M - EtEwald, 'fro') / max(norm(EtEwald, 'fro'), eps);
cosAngle = local_field_cosine(EtP3M, EtEwald);
normRatio = norm(EtP3M, 'fro') / max(norm(EtEwald, 'fro'), eps);

fprintf('\nP3M dipole cache vs Ewald reciprocal apply:\n');
fprintf('  ||Erecip Ewald||_F = %.12e\n', norm(EtEwald, 'fro'));
fprintf('  ||Erecip P3M||_F   = %.12e\n', norm(EtP3M, 'fro'));
fprintf('  norm ratio         = %.6e\n', normRatio);
fprintf('  rel err            = %.6e\n', relErr);
fprintf('  cosine             = %.6f\n', cosAngle);
fprintf('  full Ewald nK      = %d\n', opFull.info.nK);
fprintf('  no-k Ewald nK      = %d\n', opNoK.info.nK);
fprintf('  p3m mesh           = [%d %d %d]\n', cache.mesh_size);
fprintf('  p3m nK             = %d\n', cache.nK);

% Loose but meaningful early validation.
assert(cosAngle > 0.65, ...
    'P3M reciprocal dipole field should be broadly aligned with Ewald reciprocal field.');

assert(normRatio > 0.20 && normRatio < 5.0, ...
    'P3M reciprocal dipole field norm is grossly inconsistent with Ewald reciprocal field.');

assert(relErr < 1.5, ...
    'P3M reciprocal dipole field is too far from Ewald reciprocal field for smoke validation.');

io.assert_atomic_units(sys);
end

function c = local_field_cosine(A, B)
a = A(:);
b = B(:);

den = max(norm(a) * norm(b), eps);
c = dot(a, b) / den;
end

function sys = local_make_periodic_polsys()
sys = struct();

sys.site_pos = [
     3.0   2.0   2.0
    10.0   6.0   5.0
    17.0  11.0   9.0
     8.0  19.0  14.0
     6.0  15.0  12.0
];

sys.site_charge = zeros(5, 1);

sys.site_alpha = [
    0.08
    0.07
    0.06
    0.05
    0.00
];

sys.site_is_polarizable = [
    true
    true
    true
    true
    false
];

sys.n_sites = 5;
sys.thole_a = 0.39;

sys.site_type = {'X'; 'X'; 'X'; 'X'; 'X'};
sys.site_class = {'pol'; 'pol'; 'pol'; 'pol'; 'spectator'};
sys.site_label = {'p1'; 'p2'; 'p3'; 'p4'; 'n'};
sys.site_mol_id = [1; 2; 3; 4; 5];
sys.site_is_active = sys.site_is_polarizable;

sys.units = struct();
sys.units.length = 'bohr';
sys.units.alpha = 'atomic_unit';
sys.units.charge = 'elementary_charge';

sys.is_periodic = true;
sys.periodic_mode = 'periodic';

% Row-lattice convention: cart = frac * H.
sys.lattice = [
    30.0   0.0   0.0
     2.0  29.0   0.0
     1.0   3.0  28.0
];

sys.super_lattice = sys.lattice;
end