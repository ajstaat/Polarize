function [Erecip, parts] = apply_dipole_cache(cache, mu)
%APPLY_DIPOLE_CACHE Apply cached P3M reciprocal dipole field.
%
% Inputs:
%   cache  from p3m.build_dipole_cache
%   mu     N x 3 dipoles in atomic units
%
% Output:
%   Erecip N x 3 field, populated on cache.target_sites.

    if ~isstruct(cache) || ~isfield(cache, 'kind') || ...
       ~strcmp(cache.kind, 'p3m_dipole_cache')
        error('p3m:apply_dipole_cache:BadCache', ...
            'cache must come from p3m.build_dipole_cache.');
    end

    if size(mu, 1) ~= cache.nSites || size(mu, 2) ~= 3
        error('p3m:apply_dipole_cache:BadMu', ...
            'mu must be N x 3 with N matching cache.nSites.');
    end

    meshSize = cache.mesh_size;

    sourceSites = cache.source_sites;
    targetSites = cache.target_sites;

    muSource = mu(sourceSites, :);

    tAssign = tic;
    [Px, Py, Pz] = local_scatter_dipoles_cached( ...
        muSource, cache.source_stencil, meshSize);
    timeAssign = toc(tAssign);

    tFFT = tic;

    Pxk = fftn(Px);
    Pyk = fftn(Py);
    Pzk = fftn(Pz);

    PdotK = cache.kx .* Pxk + cache.ky .* Pyk + cache.kz .* Pzk;

    Exk = zeros(meshSize);
    Eyk = zeros(meshSize);
    Ezk = zeros(meshSize);

    mask = cache.mask;

    Exk(mask) = cache.Ngrid .* cache.kx(mask) .* cache.kernel(mask) .* PdotK(mask);
    Eyk(mask) = cache.Ngrid .* cache.ky(mask) .* cache.kernel(mask) .* PdotK(mask);
    Ezk(mask) = cache.Ngrid .* cache.kz(mask) .* cache.kernel(mask) .* PdotK(mask);

    Ex = real(ifftn(Exk));
    Ey = real(ifftn(Eyk));
    Ez = real(ifftn(Ezk));

    timeFFT = toc(tFFT);

    tInterp = tic;
    Etarget = local_gather_field_cached(Ex, Ey, Ez, cache.target_stencil);
    timeInterp = toc(tInterp);

    Erecip = zeros(cache.nSites, 3);
    Erecip(targetSites, :) = Etarget;

    parts = struct();
    parts.kind = 'p3m_dipole_apply';
    parts.nK = cache.nK;
    parts.nK_total_nonzero = cache.nK_total_nonzero;
    parts.mesh_size = cache.mesh_size;
    parts.assignment_order = cache.assignment_order;
    parts.alpha = cache.alpha;
    parts.kcut = cache.kcut;
    parts.use_kcut_mask = cache.use_kcut_mask;
    parts.deconvolve_assignment = cache.deconvolve_assignment;
    parts.deconvolution_floor = cache.deconvolution_floor;
    parts.lattice_convention = cache.lattice_convention;

    parts.time_assign = timeAssign;
    parts.time_fft = timeFFT;
    parts.time_interp = timeInterp;
    parts.time_total = timeAssign + timeFFT + timeInterp;
end

%% ========================================================================
% Local helpers
%% ========================================================================

function [Px, Py, Pz] = local_scatter_dipoles_cached(muSource, st, meshSize)

    Px = zeros(meshSize);
    Py = zeros(meshSize);
    Pz = zeros(meshSize);

    idx = st.linear_idx;
    w = st.weight;

    for a = 1:st.nStencil
        ia = idx(:, a);
        wa = w(:, a);

        Px = local_accum_into_grid(Px, ia, wa .* muSource(:,1));
        Py = local_accum_into_grid(Py, ia, wa .* muSource(:,2));
        Pz = local_accum_into_grid(Pz, ia, wa .* muSource(:,3));
    end
end

function grid = local_accum_into_grid(grid, idx, vals)
% Use accumarray on flattened grid. This is usually faster/cleaner than
% looping over sites with repeated indexed additions.

    flat = grid(:);
    flat = flat + accumarray(idx, vals, [numel(flat), 1], @sum, 0);
    grid = reshape(flat, size(grid));
end

function Etarget = local_gather_field_cached(Ex, Ey, Ez, st)

    Exv = Ex(:);
    Eyv = Ey(:);
    Ezv = Ez(:);

    idx = st.linear_idx;
    w = st.weight;

    n = st.nSites;

    Etarget = zeros(n, 3);

    for a = 1:st.nStencil
        ia = idx(:, a);
        wa = w(:, a);

        Etarget(:,1) = Etarget(:,1) + wa .* Exv(ia);
        Etarget(:,2) = Etarget(:,2) + wa .* Eyv(ia);
        Etarget(:,3) = Etarget(:,3) + wa .* Ezv(ia);
    end
end