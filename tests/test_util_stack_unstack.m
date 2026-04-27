function test_util_stack_unstack()
%TEST_UTIL_STACK_UNSTACK Verify stack/unstack vector conventions.
%
% This test checks the basic N x 3 <-> 3N x 1 vector convention used for
% dipoles, fields, and active-space vectors.

rng(1);

mu = randn(5, 3);

v = util.stack_xyz(mu);
mu2 = util.unstack_xyz(v);

assert(isequal(size(v), [15, 1]), ...
    'stack_xyz should produce a 3N x 1 column vector.');

assert(isequal(size(mu2), [5, 3]), ...
    'unstack_xyz should recover an N x 3 array.');

assert(norm(mu(:) - mu2(:)) < 1e-12, ...
    'unstack_xyz(stack_xyz(mu)) did not round-trip.');

% Also check a single-site edge case.
mu1 = randn(1, 3);

v1 = util.stack_xyz(mu1);
mu1b = util.unstack_xyz(v1);

assert(isequal(size(v1), [3, 1]), ...
    'Single-site stack_xyz should produce a 3 x 1 vector.');

assert(isequal(size(mu1b), [1, 3]), ...
    'Single-site unstack_xyz should recover a 1 x 3 array.');

assert(norm(mu1(:) - mu1b(:)) < 1e-12, ...
    'Single-site stack/unstack did not round-trip.');

end