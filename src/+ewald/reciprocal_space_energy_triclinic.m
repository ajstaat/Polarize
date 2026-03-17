function [Urecip, meta] = reciprocal_space_energy_triclinic(mu, r, H, alpha, kcut)
%RECIPROCAL_SPACE_ENERGY_TRICLINIC Reciprocal-space part of triclinic dipole Ewald energy.
%
% Inputs
%   mu     : N x 3 dipole vectors
%   r      : N x 3 Cartesian positions
%   H      : 3 x 3 direct lattice matrix, columns are lattice vectors
%   alpha  : Ewald screening parameter
%   kcut   : reciprocal-space cutoff in |k|
%
% Outputs
%   Urecip : reciprocal-space contribution to dipole Ewald energy
%   meta   : struct with bookkeeping info
%            .num_kvec = number of reciprocal vectors included
%
% Notes
%   - This is a direct refactor of the reference triclinic implementation.
%   - Electrostatic prefactor is 1.

if size(mu,2) ~= 3 || size(r,2) ~= 3
    error('mu and r must both be N x 3 arrays.');
end
if size(mu,1) ~= size(r,1)
    error('mu and r must have the same number of rows.');
end
if ~isequal(size(H), [3 3])
    error('H must be a 3x3 cell matrix with lattice vectors as columns.');
end
if alpha <= 0 || kcut <= 0
    error('alpha and kcut must be positive.');
end

N = size(r,1);
V = abs(det(H));

if V <= 0
    error('Cell matrix must have positive nonzero volume.');
end

G = ewald.reciprocal_lattice(H);
g1 = G(:,1);
g2 = G(:,2);
g3 = G(:,3);

% Conservative reciprocal bounds, matching your reference logic
hkx_max = ceil(kcut / norm(g1)) + 1;
hky_max = ceil(kcut / norm(g2)) + 1;
hkz_max = ceil(kcut / norm(g3)) + 1;

kvecs = zeros((2*hkx_max+1)*(2*hky_max+1)*(2*hkz_max+1), 3);
nk = 0;

for hx = -hkx_max:hkx_max
    for hy = -hky_max:hky_max
        for hz = -hkz_max:hkz_max
            if hx == 0 && hy == 0 && hz == 0
                continue;
            end

            mvec = [hx; hy; hz];
            kvec = (G * mvec).';

            if norm(kvec) <= kcut
                nk = nk + 1;
                kvecs(nk,:) = kvec;
            end
        end
    end
end

kvecs = kvecs(1:nk,:);

Urecip = 0.0;

for m = 1:nk
    kvec = kvecs(m,:);
    k2 = dot(kvec, kvec);

    pref = (2*pi / V) * exp(-k2 / (4*alpha^2)) / k2;

    Smu = 0.0 + 0.0i;
    for j = 1:N
        phase = dot(kvec, r(j,:));
        Smu = Smu + dot(kvec, mu(j,:)) * exp(1i * phase);
    end

    Urecip = Urecip + pref * abs(Smu)^2;
end

meta = struct();
meta.num_kvec = nk;

end