function [kvecs, meta] = enumerate_kvecs_triclinic(H, kcut)
%ENUMERATE_KVECS_TRICLINIC Generate reciprocal vectors with |k| <= kcut.
%
% Inputs
%   H      3x3 direct lattice matrix, columns are lattice vectors
%   kcut   reciprocal-space cutoff in |k|
%
% Outputs
%   kvecs  Nk x 3 array of reciprocal vectors
%   meta   struct with bookkeeping fields

if ~isequal(size(H), [3 3])
    error('H must be 3x3.');
end
if kcut <= 0
    error('kcut must be positive.');
end

G = ewald.reciprocal_lattice(H);
g1 = G(:,1);
g2 = G(:,2);
g3 = G(:,3);

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

meta = struct();
meta.num_kvec = nk;
meta.hkmax = [hkx_max, hky_max, hkz_max];
end