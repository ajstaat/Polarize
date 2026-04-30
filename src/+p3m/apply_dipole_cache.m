function [E, parts] = apply_dipole_cache(cache, mu)
%APPLY_DIPOLE_CACHE Apply cached P3M reciprocal dipole field.
%
% [E, parts] = p3m.apply_dipole_cache(cache, mu)
%
% Inputs
%   cache   struct from p3m.build_dipole_cache
%
%   mu      one of:
%             - Nsites x 3 full-system Cartesian dipoles
%             - Nsources x 3 source-local Cartesian dipoles
%             - 3*Nsources x 1 source-local stacked xyz vector
%
% Output
%   E       Nsites x 3 field, nonzero only on cache.target_sites
%
%   parts   diagnostic struct with mesh fields and source/target data
%
% Important
% ---------
% This applies the RECIPROCAL P3M dipole field only.
%
% It does not include:
%   - real-space dipole interactions
%   - Thole short-range correction
%   - analytic self term
%   - surface term

local_validate_cache(cache);

muSource = local_extract_source_mu(cache, mu);

M = cache.mesh_size;
order = cache.assignment_order;

Px = zeros(M);
Py = zeros(M);
Pz = zeros(M);

src = cache.source_stencil;

for a = 1:cache.nSources
    mux = muSource(a, 1);
    muy = muSource(a, 2);
    muz = muSource(a, 3);

    if mux == 0 && muy == 0 && muz == 0
        continue;
    end

    for aa = 1:order
        ii = src.i1(a, aa);
        wx = src.w1(a, aa);

        for bb = 1:order
            jj = src.i2(a, bb);
            wxy = wx * src.w2(a, bb);

            for cc = 1:order
                kk = src.i3(a, cc);
                w = wxy * src.w3(a, cc);

                Px(ii,jj,kk) = Px(ii,jj,kk) + mux * w;
                Py(ii,jj,kk) = Py(ii,jj,kk) + muy * w;
                Pz(ii,jj,kk) = Pz(ii,jj,kk) + muz * w;
            end
        end
    end
end

Pkx = fftn(Px);
Pky = fftn(Py);
Pkz = fftn(Pz);

D = cache.D;
mask = cache.maskK;

S = D.x .* Pkx + D.y .* Pky + D.z .* Pkz;

Exk = zeros(M);
Eyk = zeros(M);
Ezk = zeros(M);

% Match p3m.solve_dipole_field_spectral:
%
%   E_rec(k) = + G(k) D(k) [D(k) . P(k)]
%
% with MATLAB FFT normalization supplied by Ngrid before ifftn.
Exk(mask) = -cache.Ngrid .* D.x(mask) .* cache.influence(mask) .* S(mask);
Eyk(mask) = -cache.Ngrid .* D.y(mask) .* cache.influence(mask) .* S(mask);
Ezk(mask) = -cache.Ngrid .* D.z(mask) .* cache.influence(mask) .* S(mask);

Ex = real(ifftn(Exk));
Ey = real(ifftn(Eyk));
Ez = real(ifftn(Ezk));

Et = zeros(cache.nTargets, 3);

trg = cache.target_stencil;

for a = 1:cache.nTargets
    e = [0.0, 0.0, 0.0];

    for aa = 1:order
        ii = trg.i1(a, aa);
        wx = trg.w1(a, aa);

        for bb = 1:order
            jj = trg.i2(a, bb);
            wxy = wx * trg.w2(a, bb);

            for cc = 1:order
                kk = trg.i3(a, cc);
                w = wxy * trg.w3(a, cc);

                e(1) = e(1) + w * Ex(ii,jj,kk);
                e(2) = e(2) + w * Ey(ii,jj,kk);
                e(3) = e(3) + w * Ez(ii,jj,kk);
            end
        end
    end

    Et(a, :) = e;
end

E = zeros(cache.nSites, 3);
E(cache.target_sites, :) = Et;

parts = struct();

parts.recip = E;
parts.total = E;

parts.muSource = muSource;
parts.Mmu = sum(muSource, 1);

parts.P_total = [sum(Px(:)), sum(Py(:)), sum(Pz(:))];

parts.grid_dipoles = struct();
parts.grid_dipoles.Px = Px;
parts.grid_dipoles.Py = Py;
parts.grid_dipoles.Pz = Pz;

parts.grid_field = struct();
parts.grid_field.Ex = Ex;
parts.grid_field.Ey = Ey;
parts.grid_field.Ez = Ez;

parts.mesh_size = cache.mesh_size;
parts.assignment_order = cache.assignment_order;

parts.alpha = cache.alpha;
parts.derivative_mode = cache.derivative_mode;
parts.fd_stencil = cache.fd_stencil;
parts.influence_mode = cache.influence_mode;
parts.deconvolve_assignment = cache.deconvolve_assignment;
parts.deconvolution_floor = cache.deconvolution_floor;
parts.alias_range = cache.alias_range;

parts.source_sites = cache.source_sites;
parts.target_sites = cache.target_sites;
parts.source_mask = cache.source_mask;
parts.target_mask = cache.target_mask;

parts.nK = cache.nK;
parts.Ngrid = cache.Ngrid;

parts.lattice = cache.lattice;
parts.H = cache.H;
parts.G = cache.G;
parts.volume = cache.volume;
parts.lattice_convention = cache.lattice_convention;

parts.norm_Erecip = norm(E, 'fro');
parts.norm_muSource = norm(muSource, 'fro');
end

% =========================================================================
% Helpers
% =========================================================================

function muSource = local_extract_source_mu(cache, mu)
if isnumeric(mu) && isvector(mu)
    mu = mu(:);

    if numel(mu) ~= 3 * cache.nSources
        error('p3m:apply_dipole_cache:BadMuVectorSize', ...
            'Stacked mu vector must have length 3*cache.nSources.');
    end

    muSource = util.unstack_xyz(mu);
    return;
end

validateattributes(mu, {'numeric'}, ...
    {'2d','ncols',3,'real','finite'}, ...
    mfilename, 'mu');

if size(mu, 1) == cache.nSites
    muSource = mu(cache.source_sites, :);
elseif size(mu, 1) == cache.nSources
    muSource = mu;
else
    error('p3m:apply_dipole_cache:BadMuMatrixSize', ...
        ['mu must be either Nsites x 3, Nsources x 3, ', ...
         'or a stacked 3*Nsources vector.']);
end
end

function local_validate_cache(cache)
required = {
    'mode'
    'nSites'
    'nSources'
    'nTargets'
    'source_sites'
    'target_sites'
    'source_mask'
    'target_mask'
    'mesh_size'
    'assignment_order'
    'Ngrid'
    'source_stencil'
    'target_stencil'
    'D'
    'influence'
    'maskK'
    'alpha'
    'nK'
    'H'
    'G'
    'volume'
    'lattice'
    'lattice_convention'
};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(cache, name) || isempty(cache.(name))
        error('p3m:apply_dipole_cache:BadCache', ...
            'cache.%s is required and missing/empty.', name);
    end
end

if ~strcmp(cache.mode, 'p3m_dipole_reciprocal_cache')
    error('p3m:apply_dipole_cache:BadCacheMode', ...
        'Expected cache.mode = p3m_dipole_reciprocal_cache.');
end

if ~strcmp(cache.lattice_convention, 'project_row_H_column_G_HG_2piI')
    error('p3m:apply_dipole_cache:BadLatticeConvention', ...
        'Unexpected cache.lattice_convention "%s".', cache.lattice_convention);
end

if norm(cache.H * cache.G - 2*pi*eye(3), 'fro') > 1e-9
    error('p3m:apply_dipole_cache:BadReciprocalLattice', ...
        'Expected cache.H * cache.G = 2*pi*I.');
end

if ~isequal(size(cache.influence), cache.mesh_size)
    error('p3m:apply_dipole_cache:BadInfluenceSize', ...
        'cache.influence size must match cache.mesh_size.');
end

if ~isequal(size(cache.maskK), cache.mesh_size)
    error('p3m:apply_dipole_cache:BadMaskSize', ...
        'cache.maskK size must match cache.mesh_size.');
end
end