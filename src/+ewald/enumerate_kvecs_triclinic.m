function [kvecs, meta] = enumerate_kvecs_triclinic(H, kcut)
%ENUMERATE_KVECS_TRICLINIC Generate unique reciprocal vectors with |k| <= kcut.
%
% [kvecs, meta] = ewald.enumerate_kvecs_triclinic(H, kcut)
%
% Inputs
%   H      3x3 direct lattice matrix, columns are lattice vectors
%   kcut   reciprocal-space cutoff in |k|
%
% Outputs
%   kvecs  Nk x 3 array of unique reciprocal vectors (one from each ±k pair)
%   meta   struct with bookkeeping fields:
%       .num_kvec
%       .hkmax
%       .k2
%       .knorm
%
% Notes
%   - k = 0 is excluded.
%   - Returned k-vectors are sorted by increasing |k|.
%   - Only one representative from each ±k pair is kept.
%   - This is the correct companion to reciprocal formulas that already use
%     a factor of 2 * cos(k·r).

    validateattributes(H, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
        mfilename, 'H', 1);
    validateattributes(kcut, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'kcut', 2);

    G = ewald.reciprocal_lattice(H);

    g1 = G(:,1);
    g2 = G(:,2);
    g3 = G(:,3);

    hkx_max = ceil(kcut / norm(g1)) + 1;
    hky_max = ceil(kcut / norm(g2)) + 1;
    hkz_max = ceil(kcut / norm(g3)) + 1;

    nAlloc = (2 * hkx_max + 1) * (2 * hky_max + 1) * (2 * hkz_max + 1) - 1;
    kvecs = zeros(nAlloc, 3);
    k2 = zeros(nAlloc, 1);

    nk = 0;

    for hx = -hkx_max:hkx_max
        for hy = -hky_max:hky_max
            for hz = -hkz_max:hkz_max
                if hx == 0 && hy == 0 && hz == 0
                    continue;
                end

                % Keep only one representative of each ±k pair.
                if ~(hx > 0 || (hx == 0 && hy > 0) || (hx == 0 && hy == 0 && hz > 0))
                    continue;
                end

                mvec = [hx; hy; hz];
                kvec = (G * mvec).';
                k2val = dot(kvec, kvec);

                if k2val <= kcut^2
                    nk = nk + 1;
                    kvecs(nk, :) = kvec;
                    k2(nk) = k2val;
                end
            end
        end
    end

    kvecs = kvecs(1:nk, :);
    k2 = k2(1:nk);
    knorm = sqrt(k2);

    [knorm, perm] = sort(knorm, 'ascend');
    kvecs = kvecs(perm, :);
    k2 = k2(perm);

    meta = struct();
    meta.num_kvec = nk;
    meta.hkmax = [hkx_max, hky_max, hkz_max];
    meta.k2 = k2;
    meta.knorm = knorm;
end