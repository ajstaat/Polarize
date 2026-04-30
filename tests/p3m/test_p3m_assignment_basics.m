function test_p3m_assignment_basics()
%TEST_P3M_ASSIGNMENT_BASICS Low-level P3M assignment/interpolation tests.
%
% Covers:
%   p3m.bspline_weights_1d
%   p3m.assign_charges_bsplines
%   p3m.interpolate_field_bsplines
%
% Checks:
%   - B-spline weights are periodic, nonnegative, and sum to one
%   - charge assignment conserves total charge
%   - integer-cell translated fractional coordinates assign identically
%   - constant mesh fields interpolate exactly
%
% Dipole assignment is intentionally tested through the production cache
% path in:
%
%   tests/p3m/test_p3m_dipole_cache_apply_basics.m

rng(21);

meshSize = [8 7 6];
order = 4;

%% ------------------------------------------------------------------------
% B-spline weights
% -------------------------------------------------------------------------

M = 8;

testU = [-17.25, -1.0, -0.25, 0.0, 0.125, 1.2, 4.7, 7.99, 8.0, 19.33];

for k = 1:numel(testU)
    u = testU(k);

    [idx, w] = p3m.bspline_weights_1d(u, M, order);
    [idxShift, wShift] = p3m.bspline_weights_1d(u + 3*M, M, order);

    assert(numel(idx) == order, ...
        'B-spline index list should have length order.');
    assert(numel(w) == order, ...
        'B-spline weight list should have length order.');

    assert(all(idx >= 1 & idx <= M), ...
        'B-spline indices should be one-based periodic indices in 1:M.');

    assert(all(w >= -1e-14), ...
        'B-spline weights should be nonnegative up to roundoff.');

    assert(abs(sum(w) - 1) < 1e-13, ...
        'B-spline weights should sum to one.');

    assert(isequal(idx, idxShift), ...
        'B-spline periodic index wrapping failed for u and u+integer*M.');

    assert(norm(w - wShift, inf) < 1e-13, ...
        'B-spline periodic weights failed for u and u+integer*M.');
end

%% ------------------------------------------------------------------------
% Charge assignment conservation and periodic wrapping
% -------------------------------------------------------------------------

fracPos = [
     0.10   0.20   0.30
     0.95   0.05   0.40
    -0.20   1.15   0.75
     1.30  -0.25   2.10
];

q = [
    +0.7
    -0.2
    +0.1
    -0.6
];

rho = p3m.assign_charges_bsplines(fracPos, q, meshSize, order);

assert(isequal(size(rho), meshSize), ...
    'Charge mesh size should match meshSize.');

assert(abs(sum(rho(:)) - sum(q)) < 1e-13, ...
    'Charge assignment should conserve total charge.');

assert(all(isfinite(rho(:))), ...
    'Assigned charge mesh should be finite.');

% Translating by integer cells should not change the assigned mesh.
fracShift = fracPos + [
     1  0 -2
    -3  2  0
     0 -1  4
     2  3 -1
];

rhoShift = p3m.assign_charges_bsplines(fracShift, q, meshSize, order);

assert(norm(rho(:) - rhoShift(:), inf) < 1e-13, ...
    'Charge assignment should be invariant to integer-cell translations.');

% Empty assignment should produce a zero mesh of the requested size.
rhoEmpty = p3m.assign_charges_bsplines(zeros(0,3), zeros(0,1), meshSize, order);

assert(isequal(size(rhoEmpty), meshSize), ...
    'Empty charge assignment should return a mesh with requested size.');
assert(norm(rhoEmpty(:)) == 0, ...
    'Empty charge assignment should return all zeros.');

%% ------------------------------------------------------------------------
% Constant-field interpolation
% -------------------------------------------------------------------------

Ex = 1.25 * ones(meshSize);
Ey = -0.75 * ones(meshSize);
Ez = 0.50 * ones(meshSize);

E = p3m.interpolate_field_bsplines(fracPos, Ex, Ey, Ez, order);

Eref = repmat([1.25, -0.75, 0.50], size(fracPos, 1), 1);

assert(norm(E - Eref, inf) < 1e-13, ...
    'Constant mesh field should interpolate exactly.');

Eshift = p3m.interpolate_field_bsplines(fracShift, Ex, Ey, Ez, order);

assert(norm(Eshift - Eref, inf) < 1e-13, ...
    'Constant mesh field interpolation should be invariant to integer-cell shifts.');

Eempty = p3m.interpolate_field_bsplines(zeros(0,3), Ex, Ey, Ez, order);

assert(isequal(size(Eempty), [0 3]), ...
    'Empty field interpolation should return a 0 x 3 array.');
end