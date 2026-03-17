function [T, parts, meta] = assemble_periodic_interaction_matrix(sys, ewaldParams, scfParams)
%ASSEMBLE_PERIODIC_INTERACTION_MATRIX Assemble triclinic periodic dipole operator.
%
% Inputs
%   sys         struct with fields:
%               .site_pos
%               .super_lattice   OR .lattice
%               .site_is_polarizable
%
%   ewaldParams struct with fields:
%               .alpha
%               .rcut
%               .kcut
%               .boundary   optional, default 'tinfoil'
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

if nargin < 3
    scfParams = struct(); %#ok<NASGU>
end

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
    error('sys.site_pos is missing or empty.');
end
if ~isfield(sys, 'site_is_polarizable') || isempty(sys.site_is_polarizable)
    error('sys.site_is_polarizable is missing or empty.');
end

if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
    H = sys.super_lattice.';
elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
    % Your geom layer stores lattice vectors as rows, but the triclinic Ewald
    % code uses columns. So transpose here.
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

nSites = sys.n_sites;
polMask = logical(sys.site_is_polarizable(:));

Treal = zeros(3*nSites, 3*nSites);
Trecip = zeros(3*nSites, 3*nSites);
Tself = zeros(3*nSites, 3*nSites);
Tsurf = zeros(3*nSites, 3*nSites);

% Precompute k-vectors once
[kvecs, kmeta] = ewald.enumerate_kvecs_triclinic(H, kcut);

% Real-space image bounds
a = H(:,1); b = H(:,2); c = H(:,3);
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

        % Surface contribution applies to all pairs for vacuum
        Tsurf(ii, jj) = ewald.surface_tensor_block_dipole(H, boundary);

        % Reciprocal block
        Trecip(ii, jj) = ewald.reciprocal_space_tensor_block_triclinic( ...
            ri, rj, H, alpha, kcut, kvecs);

        % Real-space block: sum over images
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

                    Tij_real = Tij_real + ewald.real_space_tensor_block_triclinic(xvec, alpha);
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
end