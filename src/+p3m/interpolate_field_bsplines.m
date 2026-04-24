function E = interpolate_field_bsplines(Ex, Ey, Ez, fracPos, meshSize, order)
%INTERPOLATE_FIELD_BSPLINES Interpolate periodic mesh field to positions.
%
% Inputs
%   Ex,Ey,Ez  M1 x M2 x M3 mesh field components
%   fracPos   N x 3 fractional coordinates
%   meshSize  [M1 M2 M3]
%   order     B-spline interpolation order
%
% Output
%   E         N x 3 interpolated Cartesian field

    n = size(fracPos, 1);
    M = double(meshSize(:).');

    E = zeros(n, 3);

    for a = 1:n
        u = mod(fracPos(a, :), 1) .* M;

        [i1, w1] = p3m.bspline_weights_1d(u(1), M(1), order);
        [i2, w2] = p3m.bspline_weights_1d(u(2), M(2), order);
        [i3, w3] = p3m.bspline_weights_1d(u(3), M(3), order);

        ea = [0.0, 0.0, 0.0];

        for aa = 1:order
            for bb = 1:order
                for cc = 1:order
                    w = w1(aa) * w2(bb) * w3(cc);
                    ii = i1(aa);
                    jj = i2(bb);
                    kk = i3(cc);

                    ea(1) = ea(1) + w * Ex(ii,jj,kk);
                    ea(2) = ea(2) + w * Ey(ii,jj,kk);
                    ea(3) = ea(3) + w * Ez(ii,jj,kk);
                end
            end
        end

        E(a, :) = ea;
    end
end