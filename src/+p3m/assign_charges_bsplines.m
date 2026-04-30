function rho = assign_charges_bsplines(fracPos, q, meshSize, order)
%ASSIGN_CHARGES_BSPLINES Assign point charges to a periodic scalar mesh.
%
% rho = p3m.assign_charges_bsplines(fracPos, q, meshSize, order)
%
% Inputs
%   fracPos   N x 3 fractional coordinates. Any real values are accepted;
%             coordinates are wrapped periodically with mod(fracPos, 1).
%
%   q         N x 1 point charges in elementary charge / atomic units.
%
%   meshSize  1 x 3 positive integer mesh dimensions [M1 M2 M3].
%
%   order     positive integer B-spline assignment order.
%
% Output
%   rho       M1 x M2 x M3 scalar charge mesh.
%
% Notes
% -----
% rho stores assigned charges on grid nodes, not charge density divided by
% cell volume or grid-cell volume.
%
% The assignment is conservative:
%
%   sum(rho(:)) = sum(q)
%
% up to roundoff, because bspline_weights_1d normalizes each 1D stencil.

validateattributes(fracPos, {'numeric'}, ...
    {'2d', 'ncols', 3, 'real', 'finite'}, ...
    mfilename, 'fracPos', 1);

validateattributes(q, {'numeric'}, ...
    {'vector', 'real', 'finite'}, ...
    mfilename, 'q', 2);

if numel(q) ~= size(fracPos, 1)
    error('p3m:assign_charges_bsplines:SizeMismatch', ...
        'q must have one entry per row of fracPos.');
end

M = local_validate_mesh_size(meshSize);
order = local_validate_order(order);

n = size(fracPos, 1);

rho = zeros(M);

if n == 0
    return;
end

fracWrapped = mod(double(fracPos), 1);
q = double(q(:));

for a = 1:n
    u = fracWrapped(a, :) .* M;

    [i1, w1] = p3m.bspline_weights_1d(u(1), M(1), order);
    [i2, w2] = p3m.bspline_weights_1d(u(2), M(2), order);
    [i3, w3] = p3m.bspline_weights_1d(u(3), M(3), order);

    qa = q(a);

    for aa = 1:order
        ii = i1(aa);
        wx = w1(aa);

        for bb = 1:order
            jj = i2(bb);
            wxy = wx * w2(bb);

            for cc = 1:order
                kk = i3(cc);
                w = wxy * w3(cc);

                rho(ii,jj,kk) = rho(ii,jj,kk) + qa * w;
            end
        end
    end
end
end

function M = local_validate_mesh_size(meshSize)
validateattributes(meshSize, {'numeric'}, ...
    {'vector', 'numel', 3, 'integer', 'positive', 'finite'}, ...
    mfilename, 'meshSize');

M = double(meshSize(:).');
end

function order = local_validate_order(order)
validateattributes(order, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, ...
    mfilename, 'order');

order = double(order);
end