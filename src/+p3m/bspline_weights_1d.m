function [idx, w] = bspline_weights_1d(u, M, order)
%BSPLINE_WEIGHTS_1D Periodic centered cardinal B-spline weights.
%
% Inputs
%   u      grid coordinate in [0,M), not necessarily wrapped
%   M      number of grid points along this axis
%   order  assignment order. order=1 NGP, order=2 CIC, order=3/4 etc.
%
% Outputs
%   idx    1 x order, 1-based periodic grid indices
%   w      1 x order weights, sum(w)=1

    if order < 1 || order ~= round(order)
        error('p3m:bspline_weights_1d:BadOrder', ...
            'order must be a positive integer.');
    end

    u = mod(u, M);

    if order == 1
        k = round(u);
        idx = mod(k, M) + 1;
        w = 1.0;
        return;
    end

    % Centered B-spline support is [-order/2, order/2].
    % Include order integer grid points around u.
    k0 = floor(u - order/2) + 1;
    k = k0:(k0 + order - 1);

    x = u - k;
    w = zeros(1, order);

    for a = 1:order
        w(a) = local_centered_cardinal_bspline(x(a), order);
    end

    sw = sum(w);
    if sw <= 0 || ~isfinite(sw)
        error('p3m:bspline_weights_1d:BadWeights', ...
            'B-spline weights failed for u=%g M=%d order=%d.', u, M, order);
    end
    w = w ./ sw;

    idx = mod(k, M) + 1;
end

function y = local_centered_cardinal_bspline(x, p)
%CENTERED_CARDINAL_BSPLINE Cardinal B-spline of order p centered at zero.
%
% Uses:
%
%   beta_p(x) = 1/(p-1)! sum_{j=0}^{p} (-1)^j C(p,j)
%               (x + p/2 - j)_+^{p-1}
%
% Support: [-p/2, p/2].

    if abs(x) >= p/2
        % Endpoint can be exactly nonzero only for order 1; handled above.
        y = 0.0;
        return;
    end

    deg = p - 1;
    y = 0.0;

    for j = 0:p
        t = x + p/2 - j;
        if t > 0
            y = y + (-1)^j * nchoosek(p, j) * t^deg;
        end
    end

    y = y / factorial(deg);

    % Clean tiny roundoff.
    if y < 0 && y > -1e-14
        y = 0.0;
    end
end