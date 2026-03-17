function [Ureal, meta] = real_space_energy_triclinic(mu, r, H, alpha, rcut)
%REAL_SPACE_ENERGY_TRICLINIC Real-space part of triclinic dipole Ewald energy.
%
% Inputs
%   mu     : N x 3 dipole vectors
%   r      : N x 3 Cartesian positions
%   H      : 3 x 3 direct lattice matrix, columns are lattice vectors
%   alpha  : Ewald screening parameter
%   rcut   : real-space cutoff
%
% Outputs
%   Ureal  : real-space contribution to dipole Ewald energy
%   meta   : struct with bookkeeping info
%            .real_image_bounds = [nxmax nymax nzmax]
%
% Notes
%   - This is a direct refactor of the reference triclinic implementation.
%   - Electrostatic prefactor is 1.
%   - The expression matches the energy-only formulation.

if size(mu,2) ~= 3 || size(r,2) ~= 3
    error('mu and r must both be N x 3 arrays.');
end
if size(mu,1) ~= size(r,1)
    error('mu and r must have the same number of rows.');
end
if ~isequal(size(H), [3 3])
    error('H must be a 3x3 cell matrix with lattice vectors as columns.');
end
if alpha <= 0 || rcut <= 0
    error('alpha and rcut must be positive.');
end

N = size(r,1);

a = H(:,1);
b = H(:,2);
c = H(:,3);

% Conservative image bounds, matching your reference logic
nxmax = ceil(rcut / norm(a)) + 1;
nymax = ceil(rcut / norm(b)) + 1;
nzmax = ceil(rcut / norm(c)) + 1;

Ureal = 0.0;

for i = 1:N
    mui = mu(i,:);

    for j = 1:N
        muj = mu(j,:);
        rij0 = r(i,:) - r(j,:);

        for nx = -nxmax:nxmax
            for ny = -nymax:nymax
                for nz = -nzmax:nzmax

                    if i == j && nx == 0 && ny == 0 && nz == 0
                        continue;
                    end

                    nvec = [nx; ny; nz];
                    Rimg = (H * nvec).';
                    xvec = rij0 + Rimg;
                    x = norm(xvec);

                    if x == 0 || x > rcut
                        continue;
                    end

                    erfcax = erfc(alpha * x);
                    expax2 = exp(-(alpha^2) * (x^2));

                    B = erfcax / x^3 + (2*alpha/sqrt(pi)) * expax2 / x^2;
                    C = 3*erfcax / x^5 + (2*alpha/sqrt(pi)) * ...
                        (2*alpha^2 / x^2 + 3 / x^4) * expax2;

                    mui_dot_muj = dot(mui, muj);
                    mui_dot_x   = dot(mui, xvec);
                    muj_dot_x   = dot(muj, xvec);

                    Ureal = Ureal + 0.5 * ( ...
                        B * mui_dot_muj - C * mui_dot_x * muj_dot_x );
                end
            end
        end
    end
end

meta = struct();
meta.real_image_bounds = [nxmax nymax nzmax];

end