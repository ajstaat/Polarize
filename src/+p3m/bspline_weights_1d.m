function [idx, w] = bspline_weights_1d(u, M, order)
%BSPLINE_WEIGHTS_1D Periodic cardinal B-spline assignment weights in 1D.
%
% [idx, w] = p3m.bspline_weights_1d(u, M, order)
%
% Inputs
%   u       scalar mesh coordinate in grid-index units, zero-based.
%           Example: u = frac * M.
%
%   M       positive integer number of mesh points in this direction.
%
%   order   positive integer B-spline assignment order.
%           Typical P3M/PME values are 3, 4, 5, or 6.
%
% Outputs
%   idx     order x 1 one-based periodic MATLAB mesh indices in 1:M
%   w       order x 1 assignment weights, normalized to sum(w) = 1
%
% Notes
% -----
% This helper is intentionally small and dependency-free. It evaluates the
% centered cardinal B-spline using the truncated-power formula:
%
%   M_p(x) = 1/(p-1)! sum_{k=0}^{p} (-1)^k C(p,k)
%            (x + p/2 - k)_+^{p-1}
%
% with p = order. The support is centered on u and spans order grid nodes.
%
% Periodic wrapping is applied only to returned indices. The local coordinate
% inside the B-spline stencil is evaluated before wrapping.

validateattributes(u, {'numeric'}, ...
    {'scalar', 'real', 'finite'}, ...
    mfilename, 'u', 1);

validateattributes(M, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, ...
    mfilename, 'M', 2);

validateattributes(order, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, ...
    mfilename, 'order', 3);

M = double(M);
order = double(order);

% Wrap u to the primary mesh period for numerical stability.
u = mod(double(u), M);

% Centered stencil. For even order=4, this gives four nodes around u in the
% standard PME/P3M assignment sense.
left = floor(u - order/2 + 1);
nodes = left:(left + order - 1);

x = u - nodes(:);

w = local_cardinal_bspline_centered(x, order);

% Normalize defensively. The formula should already satisfy this up to
% roundoff, but normalization prevents tiny accumulated drift in tests.
sw = sum(w);
if sw <= 0 || ~isfinite(sw)
    error('p3m:bspline_weights_1d:BadWeightSum', ...
        'B-spline weights have invalid sum %.16e.', sw);
end

w = w ./ sw;

% Convert zero-based mesh nodes to one-based periodic MATLAB indices.
idx = mod(nodes(:), M) + 1;
idx = double(idx);
end

function y = local_cardinal_bspline_centered(x, order)
%LOCAL_CARDINAL_BSPLINE_CENTERED Centered cardinal B-spline of given order.

x = double(x(:));
p = double(order);

y = zeros(size(x));

shifted = x + p/2;

for k = 0:p
    z = shifted - k;
    mask = z > 0;

    if any(mask)
        y(mask) = y(mask) + ...
            (-1)^k * nchoosek(p, k) .* z(mask).^(p - 1);
    end
end

y = y ./ factorial(p - 1);

% Roundoff can produce tiny negative values near the support boundary.
tinyNeg = (y < 0) & (y > -1e-14);
y(tinyNeg) = 0;

if any(y < -1e-12)
    error('p3m:bspline_weights_1d:NegativeWeights', ...
        'Encountered significantly negative B-spline weights.');
end
end