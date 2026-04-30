function E = interpolate_field_bsplines(fracPos, Ex, Ey, Ez, order)
%INTERPOLATE_FIELD_BSPLINES Interpolate periodic mesh field to particle sites.
%
% E = p3m.interpolate_field_bsplines(fracPos, Ex, Ey, Ez, order)
%
% Inputs
%   fracPos   N x 3 fractional coordinates. Any real values are accepted;
%             coordinates are wrapped periodically with mod(fracPos, 1).
%
%   Ex,Ey,Ez  M1 x M2 x M3 mesh field components in Cartesian coordinates.
%
%   order     positive integer B-spline interpolation order.
%
% Output
%   E         N x 3 interpolated Cartesian field values.
%
% Notes
% -----
% This uses the same B-spline stencil convention as
% p3m.assign_dipoles_bsplines and should therefore be paired with that
% assignment convention in mesh operator tests.

validateattributes(fracPos, {'numeric'}, ...
    {'2d', 'ncols', 3, 'real', 'finite'}, ...
    mfilename, 'fracPos', 1);

validateattributes(Ex, {'numeric'}, ...
    {'3d', 'real', 'finite'}, ...
    mfilename, 'Ex', 2);

validateattributes(Ey, {'numeric'}, ...
    {'3d', 'real', 'finite'}, ...
    mfilename, 'Ey', 3);

validateattributes(Ez, {'numeric'}, ...
    {'3d', 'real', 'finite'}, ...
    mfilename, 'Ez', 4);

if ~isequal(size(Ex), size(Ey), size(Ez))
    error('p3m:interpolate_field_bsplines:GridSizeMismatch', ...
        'Ex, Ey, and Ez must have identical sizes.');
end

M = double(size(Ex));
if numel(M) ~= 3
    error('p3m:interpolate_field_bsplines:BadGrid', ...
        'Ex, Ey, and Ez must be 3D arrays.');
end

validateattributes(order, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, ...
    mfilename, 'order', 5);

order = double(order);

n = size(fracPos, 1);
E = zeros(n, 3);

if n == 0
    return;
end

fracWrapped = mod(double(fracPos), 1);

for a = 1:n
    u = fracWrapped(a, :) .* M;

    [i1, w1] = p3m.bspline_weights_1d(u(1), M(1), order);
    [i2, w2] = p3m.bspline_weights_1d(u(2), M(2), order);
    [i3, w3] = p3m.bspline_weights_1d(u(3), M(3), order);

    e = [0.0, 0.0, 0.0];

    for aa = 1:order
        ii = i1(aa);
        wx = w1(aa);

        for bb = 1:order
            jj = i2(bb);
            wxy = wx * w2(bb);

            for cc = 1:order
                kk = i3(cc);
                w = wxy * w3(cc);

                e(1) = e(1) + w * Ex(ii,jj,kk);
                e(2) = e(2) + w * Ey(ii,jj,kk);
                e(3) = e(3) + w * Ez(ii,jj,kk);
            end
        end
    end

    E(a, :) = e;
end
end