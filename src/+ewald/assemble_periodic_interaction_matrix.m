function [T, parts, meta] = assemble_periodic_interaction_matrix(sys, ewaldParams, scfParams)
%ASSEMBLE_PERIODIC_INTERACTION_MATRIX Assemble triclinic periodic dipole operator.
%
% Inputs
%   sys         struct with fields:
%               .site_pos
%               .site_alpha
%               .site_is_polarizable
%               .super_lattice or .lattice
%
%   ewaldParams struct with fields:
%               .alpha
%               .rcut
%               .kcut
%               .boundary              optional, default 'tinfoil'
%               .use_thole_real_space  optional, default false
%               .thole_a               optional
%
%   scfParams   currently unused, included for API compatibility
%
% Outputs
%   T           3N x 3N interaction matrix
%   parts       struct with matrices:
%               .real
%               .recip
%               .self
%               .surf
%   meta        bookkeeping info
%
% Energy convention:
%   U = 0.5 * mu_vec' * T * mu_vec
%
% Notes
%   - Optional Thole correction is applied only to the real-space part.
%   - The correction is:
%         dT = (l3/x^3) I - (3 l5/x^5) xx^T
%     added onto the bare Ewald real-space block.

if nargin < 3
    scfParams = struct(); %#ok<NASGU>
end

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
    error('sys.site_pos is missing or empty.');
end
if ~isfield(sys, 'site_is_polarizable') || isempty(sys.site_is_polarizable)
    error('sys.site_is_polarizable is missing or empty.');
end
if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
    error('sys.site_alpha is missing or empty.');
end

if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
    H = sys.super_lattice.';
elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
    H = sys.lattice.';
else
    error('Need sys.super_lattice or sys.lattice.');
end

if ~isfield(ewaldParams, 'alpha') || ~isfield(ewaldParams, 'rcut') || ~isfield(ewaldParams, 'kcut')
    error('ewaldParams must contain alpha, rcut, and kcut.');
end

alpha = ewaldParams.alpha;
rcut = ewaldParams.rcut;
kcut = ewaldParams.kcut;

boundary = 'tinfoil';
if isfield(ewaldParams, 'boundary') && ~isempty(ewaldParams.boundary)
    boundary = ewaldParams.boundary;
end

use_thole_real_space = false;
if isfield(ewaldParams, 'use_thole_real_space')
    use_thole_real_space = ewaldParams.use_thole_real_space;
end

thole_a = [];
if isfield(ewaldParams, 'thole_a') && ~isempty(ewaldParams.thole_a)
    thole_a = ewaldParams.thole_a;
elseif isfield(sys, 'thole_a') && ~isempty(sys.thole_a)
    thole_a = sys.thole_a;
end

if use_thole_real_space && isempty(thole_a)
    error('Real-space Thole correction requested, but no thole_a was provided.');
end

nSites = sys.n_sites;
polMask = logical(sys.site_is_polarizable(:));
site_alpha = sys.site_alpha(:);

Treal = zeros(3*nSites, 3*nSites);
Trecip = zeros(3*nSites, 3*nSites);
Tself = zeros(3*nSites, 3*nSites);
Tsurf = zeros(3*nSites, 3*nSites);

% Precompute k-vectors
[kvecs, kmeta] = ewald.enumerate_kvecs_triclinic(H, kcut);

% Real-space image bounds
a = H(:,1);
b = H(:,2);
c = H(:,3);

nxmax = ceil(rcut / norm(a)) + 1;
nymax = ceil(rcut / norm(b)) + 1;
nzmax = ceil(rcut / norm(c)) + 1;

for i = 1:nSites
    if ~polMask(i)
        continue;
    end

    ii = util.block3(i);
    ri = sys.site_pos(i, :);

    % Self block
    Tself(ii, ii) = ewald.self_tensor_block_dipole(alpha);

    for j = 1:nSites
        if ~polMask(j)
            continue;
        end

        jj = util.block3(j);
        rj = sys.site_pos(j, :);

        % Surface contribution
        Tsurf(ii, jj) = ewald.surface_tensor_block_dipole(H, boundary);

        % Reciprocal contribution
        Trecip(ii, jj) = ewald.reciprocal_space_tensor_block_triclinic( ...
            ri, rj, H, alpha, kcut, kvecs);

        % Real-space contribution
        Tij_real = zeros(3,3);
        rij0 = ri - rj;

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

                    Tij_block = ewald.real_space_tensor_block_triclinic(xvec, alpha);

                    if use_thole_real_space
                        Tij_block = Tij_block + ...
                            ewald.real_space_tensor_block_triclinic_thole_correction( ...
                                xvec, site_alpha(i), site_alpha(j), thole_a);
                    end

                    Tij_real = Tij_real + Tij_block;
                end
            end
        end

        Treal(ii, jj) = Tij_real;
    end
end

T = Treal + Trecip + Tself + Tsurf;

parts = struct();
parts.real = Treal;
parts.recip = Trecip;
parts.self = Tself;
parts.surf = Tsurf;

meta = struct();
meta.real_image_bounds = [nxmax nymax nzmax];
meta.num_kvec = kmeta.num_kvec;
meta.volume = abs(det(H));
meta.boundary = boundary;
meta.kcut = kcut;
meta.rcut = rcut;
meta.alpha = alpha;
meta.use_thole_real_space = use_thole_real_space;

end