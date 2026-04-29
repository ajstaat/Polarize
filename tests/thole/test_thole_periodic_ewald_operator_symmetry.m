function test_thole_periodic_ewald_operator_symmetry()
%TEST_THOLE_PERIODIC_EWALD_OPERATOR_SYMMETRY Periodic Ewald operator sanity checks.
%
% Checks:
%   1. periodic paircache operator builds
%   2. op.apply is linear
%   3. op.apply is symmetric in the energy inner product:
%
%          a' * T*b == b' * T*a
%
%   4. paircache and rowcache operators agree on op.apply
%
% This catches sign/factor/direction mistakes in real-space, reciprocal,
% self, and surface pieces.

rng(11);

sys = local_make_periodic_polsys();
Eext = local_make_external_field(sys);

scfParams = struct();
scfParams.tol = 1e-11;
scfParams.maxIter = 200;
scfParams.mixing = 0.5;
scfParams.omega = 1.0;
scfParams.verbose = false;

problem = thole.prepare_scf_problem(sys, Eext, scfParams);

baseArgs = { ...
    'Mode', 'periodic_ewald', ...
    'UseThole', true, ...
    'Rcut', 9.0, ...
    'Alpha', 0.30, ...
    'Kcut', 1.50, ...
    'Boundary', 'tinfoil', ...
    'KspaceMode', 'full', ...
    'UseMex', false, ...
    'Profile', false, ...
    'Verbose', false};

opPair = thole.make_polarization_operator( ...
    sys, problem, baseArgs{:}, ...
    'Solver', 'gmres', ...
    'Backend', 'auto');

assert(strcmp(opPair.mode, 'periodic_ewald'), ...
    'Expected periodic_ewald operator mode.');
assert(strcmp(opPair.backend, 'periodic_paircache_apply'), ...
    'Expected periodic paircache backend.');
assert(opPair.capabilities.apply, ...
    'Paircache operator should expose apply.');
assert(~opPair.capabilities.row_update, ...
    'Paircache operator should not expose production row updates.');
assert(isfield(opPair.info, 'nK') && opPair.info.nK > 0, ...
    'Periodic operator should include reciprocal k-vectors.');

nDim = 3 * problem.nPolSites;

a = randn(nDim, 1);
b = randn(nDim, 1);
c = randn(nDim, 1);

Ta = opPair.apply(a);
Tb = opPair.apply(b);
Tc = opPair.apply(c);

assert(all(isfinite(Ta)) && all(isfinite(Tb)) && all(isfinite(Tc)), ...
    'Periodic operator apply returned non-finite values.');

% Linearity.
lambda = -0.37;
lhsLin = opPair.apply(a + lambda*b);
rhsLin = Ta + lambda*Tb;

relLin = norm(lhsLin - rhsLin) / max(1, norm(rhsLin));
assert(relLin < 1e-11, ...
    'Periodic operator apply failed linearity check: rel error %.3e.', relLin);

% Symmetry.
lhsSym = a.' * Tb;
rhsSym = b.' * Ta;

relSym = abs(lhsSym - rhsSym) / max([1, abs(lhsSym), abs(rhsSym)]);
assert(relSym < 1e-10, ...
    'Periodic Ewald operator is not symmetric: rel asym %.3e.', relSym);

% A second independent random pair, to reduce chance of accidental pass.
lhsSym2 = b.' * Tc;
rhsSym2 = c.' * Tb;

relSym2 = abs(lhsSym2 - rhsSym2) / max([1, abs(lhsSym2), abs(rhsSym2)]);
assert(relSym2 < 1e-10, ...
    'Periodic Ewald operator is not symmetric on second check: rel asym %.3e.', relSym2);

% Paircache and rowcache should agree for global op.apply.
opRow = thole.make_polarization_operator( ...
    sys, problem, baseArgs{:}, ...
    'Solver', 'sor', ...
    'Backend', 'auto');

assert(strcmp(opRow.mode, 'periodic_ewald'), ...
    'Expected periodic_ewald rowcache operator mode.');
assert(strcmp(opRow.backend, 'periodic_rowcache_apply'), ...
    'Expected periodic rowcache backend.');
assert(opRow.capabilities.row_update, ...
    'Rowcache operator should expose row-update capability.');

yPair = opPair.apply(a);
yRow = opRow.apply(a);

relPairRow = norm(yPair - yRow) / max(1, norm(yPair));
assert(relPairRow < 1e-11, ...
    'Periodic paircache and rowcache apply disagree: rel error %.3e.', relPairRow);

io.assert_atomic_units(sys);
end

function sys = local_make_periodic_polsys()
sys = struct();

% Four active polarizable sites plus one neutral spectator.
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