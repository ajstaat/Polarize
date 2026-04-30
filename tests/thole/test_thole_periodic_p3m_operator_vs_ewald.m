function test_thole_periodic_p3m_operator_vs_ewald()
%TEST_THOLE_PERIODIC_P3M_OPERATOR_VS_EWALD Compare P3M operators to Ewald.
%
% Checks:
%   - make_polarization_operator routes Mode="periodic_p3m"
%   - P3M paircache backend builds for GMRES/Jacobi-style global apply
%   - P3M rowcache backend builds for SOR-style backend split
%   - paircache and rowcache global apply agree exactly
%   - rowcache apply_row fallback agrees with global apply
%   - P3M paircache operator is linear
%   - P3M paircache operator is symmetric
%   - P3M full operator is highly aligned with periodic Ewald
%   - P3M reciprocal replacement error is small relative to Ewald reciprocal
%
% Important:
%   P3M differs from Ewald only in the reciprocal piece. Therefore the
%   meaningful quantitative comparison is:
%
%       T_p3m_full - T_ewald_full
%
%   normalized by:
%
%       T_ewald_recip = T_ewald_full - T_ewald_noK
%
%   not by T_ewald_full alone, because real/self/reciprocal cancellation can
%   make ||T_ewald_full * mu|| artificially small.

rng(51);

sys = local_make_periodic_polsys();
Eext = local_make_external_field(sys);

scfParams = struct();
scfParams.tol = 1e-11;
scfParams.maxIter = 200;
scfParams.mixing = 0.5;
scfParams.omega = 1.0;
scfParams.verbose = false;

problem = thole.prepare_scf_problem(sys, Eext, scfParams);

alpha = 0.30;
rcut = 9.0;
kcut = 2.0;
boundary = 'tinfoil';
meshSize = [24 22 20];

%% ------------------------------------------------------------------------
% Ewald reference operators
% -------------------------------------------------------------------------

opEwald = thole.make_polarization_operator(sys, problem, ...
    'Mode', 'periodic_ewald', ...
    'Solver', 'gmres', ...
    'Backend', 'auto', ...
    'UseThole', true, ...
    'Rcut', rcut, ...
    'Alpha', alpha, ...
    'Kcut', kcut, ...
    'Boundary', boundary, ...
    'KspaceMode', 'full', ...
    'UseMex', false, ...
    'Profile', false, ...
    'Verbose', false);

% No-k Ewald operator. This gives real + self + surface. Subtracting it
% from the full Ewald operator extracts the reciprocal Ewald piece.
opEwaldNoK = thole.make_polarization_operator(sys, problem, ...
    'Mode', 'periodic_ewald', ...
    'Solver', 'gmres', ...
    'Backend', 'auto', ...
    'UseThole', true, ...
    'Rcut', rcut, ...
    'Alpha', alpha, ...
    'Kcut', 1.0e-12, ...
    'Boundary', boundary, ...
    'KspaceMode', 'full', ...
    'UseMex', false, ...
    'Profile', false, ...
    'Verbose', false);

assert(strcmp(opEwald.mode, 'periodic_ewald'), ...
    'Expected periodic_ewald reference operator.');
assert(opEwald.capabilities.apply, ...
    'Periodic Ewald reference operator should support apply.');

assert(isfield(opEwald.info, 'nK') && opEwald.info.nK > 0, ...
    'Full Ewald reference should have reciprocal k-vectors.');
assert(isfield(opEwaldNoK.info, 'nK') && opEwaldNoK.info.nK == 0, ...
    'No-k Ewald reference should have zero reciprocal k-vectors.');

%% ------------------------------------------------------------------------
% P3M paircache/global-apply operator
% -------------------------------------------------------------------------

opP3MPair = thole.make_polarization_operator(sys, problem, ...
    'Mode', 'periodic_p3m', ...
    'Solver', 'gmres', ...
    'Backend', 'auto', ...
    'UseThole', true, ...
    'Softening', 0.0, ...
    'Rcut', rcut, ...
    'Alpha', alpha, ...
    'Kcut', kcut, ...
    'Boundary', boundary, ...
    'MeshSize', meshSize, ...
    'AssignmentOrder', 4, ...
    'DerivativeMode', 'spectral', ...
    'InfluenceMode', 'ewald', ...
    'DeconvolveAssignment', true, ...
    'DeconvolutionFloor', 1e-8, ...
    'AliasRange', 2, ...
    'UseMex', false, ...
    'Profile', false, ...
    'Verbose', false);

assert(strcmp(opP3MPair.mode, 'periodic_p3m'), ...
    'Expected periodic_p3m pair operator mode.');
assert(strcmp(opP3MPair.backend, 'periodic_p3m_paircache_apply'), ...
    'Expected periodic_p3m_paircache_apply backend.');
assert(opP3MPair.capabilities.apply, ...
    'P3M paircache operator should support global apply.');
assert(~opP3MPair.capabilities.dense_matrix, ...
    'P3M paircache operator should not expose a dense matrix.');
assert(~opP3MPair.capabilities.row_update, ...
    'P3M paircache operator should not claim row_update capability.');

assert(isfield(opP3MPair, 'p3m_cache'), ...
    'P3M paircache operator should retain p3m_cache.');
assert(strcmp(opP3MPair.p3m_cache.mode, 'p3m_dipole_reciprocal_cache'), ...
    'P3M paircache p3m_cache should be a dipole reciprocal cache.');

assert(isfield(opP3MPair, 'periodic_p3m_cache'), ...
    'P3M paircache operator should expose periodic_p3m_cache.');
assert(isfield(opP3MPair, 'periodic_cache'), ...
    'P3M paircache operator should expose periodic_cache compatibility field.');

assert(opP3MPair.info.nK > 0, ...
    'P3M paircache operator should include reciprocal mesh modes.');
assert(opP3MPair.info.nRealEntriesDirected >= 0, ...
    'P3M paircache operator should report real-space entries.');

%% ------------------------------------------------------------------------
% P3M rowcache/SOR-oriented operator
% -------------------------------------------------------------------------

opP3MRow = thole.make_polarization_operator(sys, problem, ...
    'Mode', 'periodic_p3m', ...
    'Solver', 'sor', ...
    'Backend', 'auto', ...
    'UseThole', true, ...
    'Softening', 0.0, ...
    'Rcut', rcut, ...
    'Alpha', alpha, ...
    'Kcut', kcut, ...
    'Boundary', boundary, ...
    'MeshSize', meshSize, ...
    'AssignmentOrder', 4, ...
    'DerivativeMode', 'spectral', ...
    'InfluenceMode', 'ewald', ...
    'DeconvolveAssignment', true, ...
    'DeconvolutionFloor', 1e-8, ...
    'AliasRange', 2, ...
    'UseMex', false, ...
    'Profile', false, ...
    'Verbose', false);

assert(strcmp(opP3MRow.mode, 'periodic_p3m'), ...
    'Expected periodic_p3m row operator mode.');
assert(strcmp(opP3MRow.backend, 'periodic_p3m_rowcache_apply'), ...
    'Expected periodic_p3m_rowcache_apply backend.');
assert(opP3MRow.capabilities.apply, ...
    'P3M rowcache operator should support global apply.');
assert(~opP3MRow.capabilities.dense_matrix, ...
    'P3M rowcache operator should not expose a dense matrix.');
assert(opP3MRow.capabilities.row_update, ...
    'P3M rowcache operator should claim row_update capability.');

assert(isfield(opP3MRow, 'apply_row'), ...
    'P3M rowcache operator should expose apply_row fallback.');
assert(isfield(opP3MRow, 'p3m_cache'), ...
    'P3M rowcache operator should retain p3m_cache.');
assert(strcmp(opP3MRow.p3m_cache.mode, 'p3m_dipole_reciprocal_cache'), ...
    'P3M rowcache p3m_cache should be a dipole reciprocal cache.');

assert(isfield(opP3MRow, 'periodic_p3m_cache'), ...
    'P3M rowcache operator should expose periodic_p3m_cache.');
assert(isfield(opP3MRow, 'periodic_cache'), ...
    'P3M rowcache operator should expose periodic_cache compatibility field.');

assert(opP3MRow.info.nK > 0, ...
    'P3M rowcache operator should include reciprocal mesh modes.');
assert(opP3MRow.info.nRealEntriesDirected >= 0, ...
    'P3M rowcache operator should report real-space entries.');

%% ------------------------------------------------------------------------
% Apply tests
% -------------------------------------------------------------------------

nDim = 3 * problem.nPolSites;

a = randn(nDim, 1);
b = randn(nDim, 1);
c = randn(nDim, 1);

TaP = opP3MPair.apply(a);
TbP = opP3MPair.apply(b);
TcP = opP3MPair.apply(c);

TaRow = opP3MRow.apply(a);
TbRow = opP3MRow.apply(b);

assert(all(isfinite(TaP)) && all(isfinite(TbP)) && all(isfinite(TcP)), ...
    'P3M paircache apply returned non-finite values.');
assert(all(isfinite(TaRow)) && all(isfinite(TbRow)), ...
    'P3M rowcache apply returned non-finite values.');

relPairRowA = norm(TaP - TaRow) / max(1, norm(TaP));
relPairRowB = norm(TbP - TbRow) / max(1, norm(TbP));

assert(relPairRowA < 1e-12, ...
    'P3M paircache and rowcache apply disagree for vector A: %.3e.', relPairRowA);
assert(relPairRowB < 1e-12, ...
    'P3M paircache and rowcache apply disagree for vector B: %.3e.', relPairRowB);

% Row fallback should extract the same active-row 3-vector as global apply.
rowLocal = 2;
rowField = opP3MRow.apply_row(rowLocal, a);
TaMat = local_unstack_xyz(TaRow);

assert(norm(rowField(:) - TaMat(rowLocal, :).') < ...
        1e-12 * max(1, norm(TaMat(rowLocal, :))), ...
    'P3M rowcache apply_row fallback disagrees with global apply.');

%% ------------------------------------------------------------------------
% Linearity
% -------------------------------------------------------------------------

lambda = -0.37;

lhsLin = opP3MPair.apply(a + lambda*b);
rhsLin = TaP + lambda*TbP;

relLin = norm(lhsLin - rhsLin) / max(1, norm(rhsLin));

assert(relLin < 1e-11, ...
    'P3M paircache operator failed linearity check: rel error %.3e.', relLin);

%% ------------------------------------------------------------------------
% Symmetry
% -------------------------------------------------------------------------

lhsSym = a.' * TbP;
rhsSym = b.' * TaP;

relAsym = abs(lhsSym - rhsSym) / max([1, abs(lhsSym), abs(rhsSym)]);

assert(relAsym < 1e-10, ...
    'P3M paircache operator is not symmetric: rel asym %.3e.', relAsym);

lhsSym2 = b.' * TcP;
rhsSym2 = c.' * TbP;

relAsym2 = abs(lhsSym2 - rhsSym2) / max([1, abs(lhsSym2), abs(rhsSym2)]);

assert(relAsym2 < 1e-10, ...
    'P3M paircache operator is not symmetric on second check: rel asym %.3e.', relAsym2);

%% ------------------------------------------------------------------------
% P3M versus Ewald
% -------------------------------------------------------------------------

TaE = opEwald.apply(a);
TbE = opEwald.apply(b);

TaNoK = opEwaldNoK.apply(a);
TbNoK = opEwaldNoK.apply(b);

TaRecipE = TaE - TaNoK;
TbRecipE = TbE - TbNoK;

diffA = TaP - TaE;
diffB = TbP - TbE;

fullRelErrA = norm(diffA) / max(norm(TaE), eps);
fullRelErrB = norm(diffB) / max(norm(TbE), eps);

recipScaledErrA = norm(diffA) / max(norm(TaRecipE), eps);
recipScaledErrB = norm(diffB) / max(norm(TbRecipE), eps);

cosFullA = local_cosine(TaP, TaE);
cosFullB = local_cosine(TbP, TbE);

cosRecipDiffA = local_cosine(diffA, TaRecipE);
cosRecipDiffB = local_cosine(diffB, TbRecipE);

% The full operator should point in essentially the same direction. This is
% robust even when full-relative error is inflated by cancellation.
assert(cosFullA > 0.999, ...
    'P3M full operator should be highly aligned with Ewald for vector A.');
assert(cosFullB > 0.999, ...
    'P3M full operator should be highly aligned with Ewald for vector B.');

% The meaningful quantitative check: P3M only approximates the reciprocal
% piece, while the real/self/surface parts are shared.
assert(recipScaledErrA < 5e-3, ...
    'P3M reciprocal-scaled error for vector A is too large: %.3e.', recipScaledErrA);
assert(recipScaledErrB < 5e-3, ...
    'P3M reciprocal-scaled error for vector B is too large: %.3e.', recipScaledErrB);

io.assert_atomic_units(sys);
end

% =========================================================================
% Helpers
% =========================================================================

function c = local_cosine(a, b)
den = max(norm(a) * norm(b), eps);
c = dot(a(:), b(:)) / den;
end

function X = local_unstack_xyz(v)
v = v(:);

if mod(numel(v), 3) ~= 0
    error('test_thole_periodic_p3m_operator_vs_ewald:BadVectorLength', ...
        'Vector length must be divisible by 3.');
end

n = numel(v) / 3;

X = zeros(n, 3);
X(:, 1) = v(1:3:end);
X(:, 2) = v(2:3:end);
X(:, 3) = v(3:3:end);
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

function Eext = local_make_external_field(sys)
Eext = zeros(sys.n_sites, 3);

Eext(1, :) = [ 1.0e-3, -0.5e-3,  0.2e-3];
Eext(2, :) = [-0.4e-3,  0.8e-3, -0.1e-3];
Eext(3, :) = [ 0.3e-3,  0.1e-3,  0.5e-3];
Eext(4, :) = [-0.2e-3, -0.3e-3,  0.4e-3];
end