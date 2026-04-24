function [kvecs, meta] = enumerate_kvecs_from_lattice(latOrSysOrH, kcut)
%ENUMERATE_KVECS_FROM_LATTICE Generate reciprocal vectors using project convention.
%
% Project lattice convention:
%
%   H has direct lattice vectors as rows:
%
%       r_cart = f_frac * H
%
%   G has reciprocal lattice vectors as columns:
%
%       H * G = 2*pi*I
%
%   For integer reciprocal index m_col:
%
%       k_col = G * m_col
%       k_row = k_col.'
%
% This function returns one representative from each +/- k pair, excluding
% k=0, with |k| <= kcut. Returned kvecs are Nk x 3 row vectors.

    validateattributes(kcut, {'double'}, ...
        {'scalar','real','finite','positive'}, ...
        mfilename, 'kcut', 2);

    lat = local_as_lattice(latOrSysOrH);

    G = lat.G;

    g1 = G(:,1);
    g2 = G(:,2);
    g3 = G(:,3);

    hkx_max = ceil(kcut / norm(g1)) + 1;
    hky_max = ceil(kcut / norm(g2)) + 1;
    hkz_max = ceil(kcut / norm(g3)) + 1;

    nAlloc = (2*hkx_max + 1) * (2*hky_max + 1) * (2*hkz_max + 1) - 1;

    kvecs = zeros(nAlloc, 3);
    k2 = zeros(nAlloc, 1);
    hkl = zeros(nAlloc, 3);

    nk = 0;

    for hx = -hkx_max:hkx_max
        for hy = -hky_max:hky_max
            for hz = -hkz_max:hkz_max
                if hx == 0 && hy == 0 && hz == 0
                    continue;
                end

                % Keep one representative of each +/- k pair.
                if ~(hx > 0 || (hx == 0 && hy > 0) || ...
                            (hx == 0 && hy == 0 && hz > 0))
                    continue;
                end

                m = [hx; hy; hz];
                krow = (G * m).';
                k2val = dot(krow, krow);

                if k2val <= kcut^2
                    nk = nk + 1;
                    kvecs(nk, :) = krow;
                    k2(nk) = k2val;
                    hkl(nk, :) = [hx hy hz];
                end
            end
        end
    end

    kvecs = kvecs(1:nk, :);
    k2 = k2(1:nk);
    hkl = hkl(1:nk, :);

    knorm = sqrt(k2);

    [knorm, perm] = sort(knorm, 'ascend');
    kvecs = kvecs(perm, :);
    k2 = k2(perm);
    hkl = hkl(perm, :);

    meta = struct();
    meta.num_kvec = nk;
    meta.hkmax = [hkx_max, hky_max, hkz_max];
    meta.k2 = k2;
    meta.knorm = knorm;
    meta.hkl = hkl;
    meta.convention = 'project_row_H_column_G_HG_2piI';
end

function lat = local_as_lattice(x)
    if isstruct(x) && isfield(x, 'H') && isfield(x, 'G')
        lat = x;
    else
        lat = geom.get_lattice(x);
    end
end