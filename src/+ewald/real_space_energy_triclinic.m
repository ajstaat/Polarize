function [Ureal, meta] = real_space_energy_triclinic(mu, r, H, alpha, rcut)
%REAL_SPACE_ENERGY_TRICLINIC Real-space part of triclinic dipole Ewald energy.
%
% [Ureal, meta] = ewald.real_space_energy_triclinic(mu, r, H, alpha, rcut)
%
% Inputs
%   mu     N x 3 dipole vectors
%   r      N x 3 Cartesian positions
%   H      3 x 3 direct lattice matrix, columns are lattice vectors
%   alpha  Ewald screening parameter
%   rcut   real-space cutoff
%
% Outputs
%   Ureal  real-space contribution to dipole Ewald energy
%   meta   struct with bookkeeping info
%          .real_image_bounds = [nxmax nymax nzmax]
%          .num_image_shifts
%
% Notes
%   - Uses i < j pair symmetry for the central cell.
%   - Includes self-image terms (i == j, R ~= 0).
%   - Electrostatic prefactor is 1.

    validateattributes(mu, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
        mfilename, 'mu', 1);
    validateattributes(r, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
        mfilename, 'r', 2);
    if size(mu, 1) ~= size(r, 1)
        error('ewald:real_space_energy_triclinic:SizeMismatch', ...
            'mu and r must have the same number of rows.');
    end
    validateattributes(H, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
        mfilename, 'H', 3);
    validateattributes(alpha, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'alpha', 4);
    validateattributes(rcut, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'rcut', 5);

    N = size(r, 1);

    a = H(:,1);
    b = H(:,2);
    c = H(:,3);

    nxmax = ceil(rcut / norm(a)) + 1;
    nymax = ceil(rcut / norm(b)) + 1;
    nzmax = ceil(rcut / norm(c)) + 1;

    % Enumerate all image shifts once.
    nAlloc = (2*nxmax + 1) * (2*nymax + 1) * (2*nzmax + 1);
    Rshifts = zeros(nAlloc, 3);
    isZeroShift = false(nAlloc, 1);

    nR = 0;
    for nx = -nxmax:nxmax
        for ny = -nymax:nymax
            for nz = -nzmax:nzmax
                nR = nR + 1;
                nvec = [nx; ny; nz];
                Rshifts(nR, :) = (H * nvec).';
                isZeroShift(nR) = (nx == 0 && ny == 0 && nz == 0);
            end
        end
    end

    Rshifts = Rshifts(1:nR, :);
    isZeroShift = isZeroShift(1:nR);

    Ureal = 0.0;
    alpha2 = alpha^2;
    twoAlphaOverSqrtPi = 2 * alpha / sqrt(pi);

    % Off-diagonal pairs: use i < j and sum over all image shifts.
    for i = 1:(N-1)
        mui = mu(i, :);
        ri = r(i, :);

        for j = (i+1):N
            muj = mu(j, :);
            rij0 = ri - r(j, :);

            xvec = rij0 + Rshifts;                % nR x 3
            x2 = sum(xvec.^2, 2);
            x = sqrt(x2);

            keep = (x > 0) & (x <= rcut);
            if ~any(keep)
                continue;
            end

            xvec = xvec(keep, :);
            x2 = x2(keep);
            x = x(keep);

            erfcax = erfc(alpha * x);
            expax2 = exp(-alpha2 * x2);

            invx2 = 1 ./ x2;
            invx = 1 ./ x;
            invx3 = invx .* invx2;
            invx5 = invx3 .* invx2;

            B = erfcax .* invx3 + twoAlphaOverSqrtPi * expax2 .* invx2;
            C = 3 * erfcax .* invx5 + twoAlphaOverSqrtPi * ...
                (2 * alpha2 * invx2 + 3 * invx5 ./ invx) .* expax2;
            % Note: invx5 ./ invx = 1/x^4
            % Equivalent to 2*a^2/x^2 + 3/x^4

            mui_dot_muj = dot(mui, muj);
            mui_dot_x = xvec * mui.';
            muj_dot_x = xvec * muj.';

            Ureal = Ureal + sum(B * mui_dot_muj - C .* mui_dot_x .* muj_dot_x);
        end
    end

    % Self-image terms: i == j, R ~= 0 only.
    nonzeroShift = ~isZeroShift;
    Rnz = Rshifts(nonzeroShift, :);

    if ~isempty(Rnz)
        x2 = sum(Rnz.^2, 2);
        x = sqrt(x2);
        keep = (x > 0) & (x <= rcut);

        if any(keep)
            Rnz = Rnz(keep, :);
            x2 = x2(keep);
            x = x(keep);

            erfcax = erfc(alpha * x);
            expax2 = exp(-alpha2 * x2);

            invx2 = 1 ./ x2;
            invx = 1 ./ x;
            invx3 = invx .* invx2;
            invx5 = invx3 .* invx2;

            B = erfcax .* invx3 + twoAlphaOverSqrtPi * expax2 .* invx2;
            C = 3 * erfcax .* invx5 + twoAlphaOverSqrtPi * ...
                (2 * alpha2 * invx2 + 3 * invx5 ./ invx) .* expax2;

            for i = 1:N
                mui = mu(i, :);
                mui2 = dot(mui, mui);
                mui_dot_x = Rnz * mui.';
                Ureal = Ureal + 0.5 * sum(B * mui2 - C .* (mui_dot_x.^2));
            end
        end
    end

    meta = struct();
    meta.real_image_bounds = [nxmax, nymax, nzmax];
    meta.num_image_shifts = nR;
end