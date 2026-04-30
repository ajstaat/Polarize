function test_p3m_kgrid_and_charge_field_smoke()
%TEST_P3M_KGRID_AND_CHARGE_FIELD_SMOKE Basic spectral charge-field tests.
%
% Covers:
%   p3m.make_kgrid
%   p3m.assign_charges_bsplines
%   p3m.solve_charge_field_spectral
%   p3m.interpolate_field_bsplines
%
% Checks:
%   - k-grid uses project row-lattice convention
%   - zero charge mesh gives zero field
%   - uniform charge mesh gives zero field because k=0 is removed
%   - neutral +q/-q source produces finite nonzero mesh/interpolated fields
%   - spectral and finite-difference derivative modes both run and return
%     finite fields
%
% This is intentionally a smoke test, not a quantitative P3M-vs-Ewald
% validation. That comes after the high-level P3M external-field driver.

rng(22);

sys = local_make_periodic_sys();
lat = geom.get_lattice(sys);

meshSize = [12 10 8];
order = 4;
alpha = 0.30;

%% ------------------------------------------------------------------------
% k-grid convention
% -------------------------------------------------------------------------

kg = p3m.make_kgrid(lat, meshSize);

assert(isequal(kg.meshSize, meshSize), ...
    'k-grid meshSize should match requested mesh.');
assert(strcmp(kg.convention, 'project_row_H_column_G_HG_2piI'), ...
    'k-grid should advertise canonical project lattice convention.');
assert(norm(kg.H * kg.G - 2*pi*eye(3), 'fro') < 1e-12, ...
    'k-grid lattice should satisfy H*G = 2*pi*I.');
assert(abs(kg.volume - lat.volume) < 1e-12 * max(1, lat.volume), ...
    'k-grid volume should match lattice volume.');

zeroIdx = (kg.m1 == 0) & (kg.m2 == 0) & (kg.m3 == 0);
assert(nnz(zeroIdx) == 1, ...
    'FFT k-grid should contain exactly one zero mode.');
assert(abs(kg.k2(zeroIdx)) < 1e-14, ...
    'Zero FFT mode should have k2 = 0.');

nonzero = ~zeroIdx;
assert(all(kg.k2(nonzero) > 0), ...
    'All nonzero FFT modes should have positive k2.');

%% ------------------------------------------------------------------------
% zero charge mesh gives zero field
% -------------------------------------------------------------------------

rho0 = zeros(meshSize);

opts = struct();
opts.alpha = alpha;
opts.assignment_order = order;
opts.derivative_mode = 'spectral';
opts.influence_mode = 'ewald';
opts.deconvolve_assignment = true;

[Ex0, Ey0, Ez0, info0] = p3m.solve_charge_field_spectral(rho0, lat, opts);

assert(norm(Ex0(:)) == 0 && norm(Ey0(:)) == 0 && norm(Ez0(:)) == 0, ...
    'Zero charge mesh should produce exactly zero spectral field.');

assert(info0.nK == prod(meshSize) - 1, ...
    'Spectral charge solve should use all nonzero FFT modes.');

%% ------------------------------------------------------------------------
% uniform charge mesh gives zero field because k=0 is removed
% -------------------------------------------------------------------------

rhoUniform = ones(meshSize);

[ExU, EyU, EzU] = p3m.solve_charge_field_spectral(rhoUniform, lat, opts);

assert(norm(ExU(:)) < 1e-12, ...
    'Uniform charge mesh should produce zero Ex after removing k=0.');
assert(norm(EyU(:)) < 1e-12, ...
    'Uniform charge mesh should produce zero Ey after removing k=0.');
assert(norm(EzU(:)) < 1e-12, ...
    'Uniform charge mesh should produce zero Ez after removing k=0.');

%% ------------------------------------------------------------------------
% neutral point charges produce finite nonzero field
% -------------------------------------------------------------------------

fracSource = [
    0.20  0.30  0.40
    0.65  0.55  0.35
];

q = [+1.0; -1.0];

rho = p3m.assign_charges_bsplines(fracSource, q, meshSize, order);

assert(abs(sum(rho(:))) < 1e-13, ...
    'Neutral assigned charge mesh should have near-zero total charge.');

[Ex, Ey, Ez, info] = p3m.solve_charge_field_spectral(rho, lat, opts);

assert(all(isfinite(Ex(:))) && all(isfinite(Ey(:))) && all(isfinite(Ez(:))), ...
    'Spectral charge field should be finite.');
assert(norm([Ex(:); Ey(:); Ez(:)]) > 0, ...
    'Neutral separated charges should produce nonzero mesh field.');

assert(strcmp(info.derivative_mode, 'spectral'), ...
    'Expected spectral derivative mode metadata.');
assert(strcmp(info.influence_mode, 'ewald'), ...
    'Expected Ewald influence mode metadata.');
assert(abs(info.total_charge_mesh - sum(q)) < 1e-13, ...
    'solve_charge_field_spectral should report conserved total mesh charge.');

fracTarget = [
    0.10  0.10  0.10
    0.40  0.50  0.60
    0.90  0.20  0.70
];

Etarg = p3m.interpolate_field_bsplines(fracTarget, Ex, Ey, Ez, order);

assert(isequal(size(Etarg), [size(fracTarget,1), 3]), ...
    'Interpolated target field should be Ntarget x 3.');
assert(all(isfinite(Etarg(:))), ...
    'Interpolated target field should be finite.');
assert(norm(Etarg, 'fro') > 0, ...
    'Interpolated target field should be nonzero for separated neutral charges.');

%% ------------------------------------------------------------------------
% finite-difference derivative mode smoke test
% -------------------------------------------------------------------------

optsFD = opts;
optsFD.derivative_mode = 'finite_difference';
optsFD.fd_stencil = 'central2';
optsFD.influence_mode = 'fd_least_squares';

[ExFD, EyFD, EzFD, infoFD] = p3m.solve_charge_field_spectral(rho, lat, optsFD);

assert(all(isfinite(ExFD(:))) && all(isfinite(EyFD(:))) && all(isfinite(EzFD(:))), ...
    'Finite-difference charge field should be finite.');
assert(norm([ExFD(:); EyFD(:); EzFD(:)]) > 0, ...
    'Finite-difference charge field should be nonzero.');

assert(strcmp(infoFD.derivative_mode, 'finite_difference'), ...
    'Expected finite_difference derivative mode metadata.');
assert(strcmp(infoFD.influence_mode, 'fd_least_squares'), ...
    'Expected fd_least_squares influence mode metadata.');

% The FD and spectral fields are not expected to be identical, but they
% should point to broadly similar physics for this smoke setup.
dotFields = dot([Ex(:); Ey(:); Ez(:)], [ExFD(:); EyFD(:); EzFD(:)]);
cosAngle = dotFields / max(norm([Ex(:); Ey(:); Ez(:)]) * norm([ExFD(:); EyFD(:); EzFD(:)]), eps);

assert(cosAngle > 0.5, ...
    'Spectral and FD charge fields should be broadly aligned in this smoke test.');
end

function sys = local_make_periodic_sys()
sys = struct();

sys.site_pos = [
     3.0   2.0   2.0
    10.0   6.0   5.0
    17.0  11.0   9.0
     8.0  19.0  14.0
];

sys.site_charge = zeros(4, 1);
sys.site_alpha = 0.05 * ones(4, 1);
sys.site_is_polarizable = true(4, 1);
sys.site_is_active = sys.site_is_polarizable;

sys.n_sites = 4;
sys.thole_a = 0.39;

sys.site_type = {'X'; 'X'; 'X'; 'X'};
sys.site_class = {'pol'; 'pol'; 'pol'; 'pol'};
sys.site_label = {'p1'; 'p2'; 'p3'; 'p4'};
sys.site_mol_id = (1:4).';

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