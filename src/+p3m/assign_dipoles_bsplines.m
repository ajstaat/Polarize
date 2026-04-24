function [Px, Py, Pz] = assign_dipoles_bsplines(fracPos, mu, meshSize, order)
%ASSIGN_DIPOLES_BSPLINES Assign point dipoles to periodic vector mesh.
%
% Inputs
%   fracPos   N x 3 fractional coordinates, any real values accepted
%   mu        N x 3 dipoles in Cartesian atomic units
%   meshSize  1 x 3 grid dimensions [M1 M2 M3]
%   order     B-spline assignment order
%
% Outputs
%   Px,Py,Pz  M1 x M2 x M3 mesh components of assigned dipole moment
%
% Notes
%   These are assigned dipole moments on grid nodes, not divided by volume.

    n = size(fracPos, 1);

    if size(fracPos, 2) ~= 3
        error('p3m:assign_dipoles_bsplines:BadFrac', ...
            'fracPos must be N x 3.');
    end
    if size(mu, 1) ~= n || size(mu, 2) ~= 3
        error('p3m:assign_dipoles_bsplines:BadMu', ...
            'mu must be N x 3.');
    end

    M = double(meshSize(:).');
    if numel(M) ~= 3 || any(M < 1) || any(M ~= round(M))
        error('p3m:assign_dipoles_bsplines:BadMesh', ...
            'meshSize must be positive integer [M1 M2 M3].');
    end

    Px = zeros(M);
    Py = zeros(M);
    Pz = zeros(M);

    for a = 1:n
        u = mod(fracPos(a, :), 1) .* M;

        [i1, w1] = p3m.bspline_weights_1d(u(1), M(1), order);
        [i2, w2] = p3m.bspline_weights_1d(u(2), M(2), order);
        [i3, w3] = p3m.bspline_weights_1d(u(3), M(3), order);

        mux = mu(a, 1);
        muy = mu(a, 2);
        muz = mu(a, 3);

        for aa = 1:order
            for bb = 1:order
                for cc = 1:order
                    w = w1(aa) * w2(bb) * w3(cc);

                    ii = i1(aa);
                    jj = i2(bb);
                    kk = i3(cc);

                    Px(ii,jj,kk) = Px(ii,jj,kk) + mux * w;
                    Py(ii,jj,kk) = Py(ii,jj,kk) + muy * w;
                    Pz(ii,jj,kk) = Pz(ii,jj,kk) + muz * w;
                end
            end
        end
    end
end