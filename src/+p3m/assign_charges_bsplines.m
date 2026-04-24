function rho = assign_charges_bsplines(fracPos, q, meshSize, order)
%ASSIGN_CHARGES_BSPLINES Assign charges to periodic mesh.
%
% Inputs
%   fracPos   N x 3 fractional coordinates, any real values accepted
%   q         N x 1 charges
%   meshSize  1 x 3 grid dimensions [M1 M2 M3]
%   order     B-spline assignment order
%
% Output
%   rho       M1 x M2 x M3 grid containing assigned charge, not density

    q = q(:);
    n = size(fracPos, 1);

    if size(fracPos, 2) ~= 3
        error('p3m:assign_charges_bsplines:BadFrac', ...
            'fracPos must be N x 3.');
    end
    if numel(q) ~= n
        error('p3m:assign_charges_bsplines:BadCharge', ...
            'q length must match fracPos rows.');
    end

    M = double(meshSize(:).');
    if numel(M) ~= 3 || any(M < 1) || any(M ~= round(M))
        error('p3m:assign_charges_bsplines:BadMesh', ...
            'meshSize must be positive integer [M1 M2 M3].');
    end

    rho = zeros(M);

    for a = 1:n
        u = mod(fracPos(a, :), 1) .* M;

        [i1, w1] = p3m.bspline_weights_1d(u(1), M(1), order);
        [i2, w2] = p3m.bspline_weights_1d(u(2), M(2), order);
        [i3, w3] = p3m.bspline_weights_1d(u(3), M(3), order);

        qa = q(a);

        for aa = 1:order
            for bb = 1:order
                for cc = 1:order
                    rho(i1(aa), i2(bb), i3(cc)) = ...
                        rho(i1(aa), i2(bb), i3(cc)) + ...
                        qa * w1(aa) * w2(bb) * w3(cc);
                end
            end
        end
    end
end