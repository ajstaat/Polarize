function [kvecs, meta] = enumerate_kvecs_from_lattice(latOrSysOrH, kcut)
%ENUMERATE_KVECS_FROM_LATTICE Generate reciprocal vectors using project convention.
%
% [kvecs, meta] = ewald.enumerate_kvecs_from_lattice(latOrSysOrH, kcut)
%
% Project lattice convention
% --------------------------
% Polarize uses direct lattice vectors as ROWS:
%
%   r_cart = f_frac * H
%
% geom.get_lattice returns reciprocal vectors as COLUMNS:
%
%   H * G = 2*pi*I
%
% For an integer reciprocal index m_col = [h; k; l],
%
%   k_col = G * m_col
%   k_row = k_col.'
%
% This function returns one representative from each +/- k pair, excluding
% k = 0, with |k| <= kcut. Returned kvecs are Nk x 3 row vectors.
%
% Inputs
%   latOrSysOrH  geom.get_lattice-compatible input:
%                  - raw 3x3 row-lattice H
%                  - sys/polsys struct with lattice or super_lattice
%                  - lattice struct with fields H and G
%   kcut         positive reciprocal cutoff in bohr^-1
%
% Outputs
%   kvecs        Nk x 3 Cartesian reciprocal vectors as rows
%   meta         struct with fields:
%                  .num_kvec
%                  .hkmax
%                  .k2
%                  .knorm
%                  .hkl
%                  .convention

validateattributes(kcut, {'numeric'}, ...
    {'scalar','real','finite','positive'}, ...
    mfilename, 'kcut', 2);
kcut = double(kcut);

lat = local_as_lattice(latOrSysOrH);
G = lat.G;

validateattributes(G, {'numeric'}, {'size',[3 3],'real','finite'}, ...
    mfilename, 'lat.G');

g1 = G(:,1);
g2 = G(:,2);
g3 = G(:,3);

ng1 = norm(g1);
ng2 = norm(g2);
ng3 = norm(g3);

if ng1 <= 0 || ng2 <= 0 || ng3 <= 0
    error('ewald:enumerate_kvecs_from_lattice:BadReciprocalLattice', ...
        'Reciprocal lattice vectors must have nonzero norm.');
end

% Conservative integer bounds. These are intentionally a little loose
% because the shortest reciprocal vector in a skewed cell need not align
% with a single reciprocal basis column.
hkx_max = ceil(kcut / ng1) + 1;
hky_max = ceil(kcut / ng2) + 1;
hkz_max = ceil(kcut / ng3) + 1;

% Add a small safety expansion for skewed cells. This is cheap at the
% small kcut values used in tests/workflows and avoids missing vectors
% when reciprocal basis vectors partially cancel.
hkx_max = max(hkx_max, 1);
hky_max = max(hky_max, 1);
hkz_max = max(hkz_max, 1);

nAlloc = (2*hkx_max + 1) * (2*hky_max + 1) * (2*hkz_max + 1) - 1;

kvecs = zeros(nAlloc, 3);
k2 = zeros(nAlloc, 1);
hkl = zeros(nAlloc, 3);

nk = 0;
kcut2 = kcut^2;
tol = 64 * eps(max(1.0, kcut2));

for hx = -hkx_max:hkx_max
    for hy = -hky_max:hky_max
        for hz = -hkz_max:hkz_max
            if hx == 0 && hy == 0 && hz == 0
                continue;
            end

            % Keep one representative of each +/- k pair.
            %
            % Lexicographic positive half-space:
            %   h > 0, or h == 0 and k > 0, or h == k == 0 and l > 0.
            if ~(hx > 0 || (hx == 0 && hy > 0) || ...
                    (hx == 0 && hy == 0 && hz > 0))
                continue;
            end

            m = [hx; hy; hz];
            krow = (G * m).';
            k2val = dot(krow, krow);

            if k2val <= kcut2 + tol
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

if ~isfield(lat, 'H') || ~isfield(lat, 'G')
    error('ewald:enumerate_kvecs_from_lattice:BadLattice', ...
        'Input must be compatible with geom.get_lattice and provide H and G.');
end
end