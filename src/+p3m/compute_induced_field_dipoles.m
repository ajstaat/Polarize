function [E, parts] = compute_induced_field_dipoles(sys, mu, opts)
%COMPUTE_INDUCED_FIELD_DIPOLES P3M reciprocal/surface field from dipoles.
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
% Computes only:
%
%   E = E_recip,P3M + E_surf
%
% It does NOT compute the real-space dipole field and does NOT include the
% Ewald self term.

    io.assert_atomic_units(sys);

    if nargin < 3
        opts = struct();
    end

    if ~isfield(opts, 'ewald') || isempty(opts.ewald)
        error('p3m:compute_induced_field_dipoles:MissingEwald', ...
            'opts.ewald is required.');
    end
    if ~isfield(opts, 'mesh_size') || isempty(opts.mesh_size)
        error('p3m:compute_induced_field_dipoles:MissingMeshSize', ...
            'opts.mesh_size is required.');
    end
    if ~isfield(opts, 'assignment_order') || isempty(opts.assignment_order)
        error('p3m:compute_induced_field_dipoles:MissingAssignmentOrder', ...
            'opts.assignment_order is required.');
    end

    alpha = local_require_struct_field(opts.ewald, 'alpha');

    boundary = 'tinfoil';
    if isfield(opts.ewald, 'boundary') && ~isempty(opts.ewald.boundary)
        boundary = lower(char(string(opts.ewald.boundary)));
    end
    if ~ismember(boundary, {'tinfoil','vacuum'})
        error('p3m:compute_induced_field_dipoles:BadBoundary', ...
            'boundary must be ''tinfoil'' or ''vacuum''.');
    end

    meshSize = double(opts.mesh_size(:).');
    if numel(meshSize) ~= 3 || any(meshSize < 1) || any(meshSize ~= round(meshSize))
        error('p3m:compute_induced_field_dipoles:BadMeshSize', ...
            'opts.mesh_size must be positive integer [M1 M2 M3].');
    end

    order = opts.assignment_order;
    validateattributes(order, {'numeric'}, ...
        {'scalar','real','finite','positive','integer'}, ...
        mfilename, 'opts.assignment_order');
    order = double(order);

    pos = sys.site_pos;
    nSites = size(pos, 1);

    if size(mu, 1) ~= nSites || size(mu, 2) ~= 3
        error('p3m:compute_induced_field_dipoles:BadMu', ...
            'mu must be N x 3.');
    end

    targetMask = local_get_mask(opts, 'target_mask', true(nSites, 1));

    if isfield(opts, 'source_mask') && ~isempty(opts.source_mask)
        sourceMask = logical(opts.source_mask(:));
    else
        sourceMask = vecnorm(mu, 2, 2) > 0;
    end

    if numel(targetMask) ~= nSites || numel(sourceMask) ~= nSites
        error('p3m:compute_induced_field_dipoles:BadMask', ...
            'target_mask and source_mask must have length N.');
    end

    lat = geom.get_lattice(sys);
    H = lat.H;
    V = lat.volume;

    if V <= 1e-14
        error('p3m:compute_induced_field_dipoles:BadVolume', ...
            'Cell volume is nonpositive.');
    end

    verbose = local_get_opt(opts, 'verbose', false);

    targetSites = find(targetMask);
    sourceSites = find(sourceMask);

    if verbose
        fprintf('p3m.compute_induced_field_dipoles:\n');
        fprintf('  lattice convention = H_rows_G_columns_HG_2piI\n');
        fprintf('  mesh_size         = [%d %d %d]\n', meshSize);
        fprintf('  order             = %d\n', order);
        fprintf('  alpha             = %.8g\n', alpha);
        fprintf('  boundary          = %s\n', boundary);
        fprintf('  nTargets          = %d\n', numel(targetSites));
        fprintf('  nSources          = %d\n', numel(sourceSites));
    end

    E = zeros(nSites, 3);
    Erecip = zeros(nSites, 3);
    Esurf = zeros(nSites, 3);

    if isempty(targetSites) || isempty(sourceSites)
        parts = local_empty_parts(sys, targetMask, sourceMask, alpha, boundary);
        return;
    end

    fracSource = pos(sourceSites, :) / H;
    fracTarget = pos(targetSites, :) / H;

    muSource = mu(sourceSites, :);

    tMesh = tic;

    [Px, Py, Pz] = p3m.assign_dipoles_bsplines( ...
        fracSource, muSource, meshSize, order);

    solveOpts = struct();
    solveOpts.alpha = alpha;
    solveOpts.assignment_order = order;
    solveOpts.deconvolve_assignment = local_get_opt(opts, 'deconvolve_assignment', true);
    solveOpts.deconvolution_floor = local_get_opt(opts, 'deconvolution_floor', 1e-8);
    solveOpts.use_kcut_mask = local_get_opt(opts, 'use_kcut_mask', false);

    if solveOpts.use_kcut_mask
        solveOpts.kcut = local_get_opt(opts, 'kcut', local_get_opt(opts.ewald, 'kcut', []));
        if isempty(solveOpts.kcut)
            error('p3m:compute_induced_field_dipoles:MissingKcut', ...
                'kcut is required when use_kcut_mask=true.');
        end
    end

    [ExGrid, EyGrid, EzGrid, solveInfo] = ...
        p3m.solve_dipole_field_spectral(Px, Py, Pz, lat, solveOpts);

    ErecipTarget = p3m.interpolate_field_bsplines( ...
        ExGrid, EyGrid, EzGrid, fracTarget, meshSize, order);

    Erecip(targetSites, :) = ErecipTarget;

    timeMesh = toc(tMesh);

    Mmu = sum(muSource, 1);

    switch boundary
        case 'tinfoil'
            Esurf_mu = [0.0, 0.0, 0.0];
            surf_coeff = 0.0;

        case 'vacuum'
            surf_coeff = 4*pi/(3*V);
            Esurf_mu = -surf_coeff * Mmu;
            Esurf(targetSites, :) = repmat(Esurf_mu, numel(targetSites), 1);

        otherwise
            error('p3m:compute_induced_field_dipoles:BadBoundaryInternal', ...
                'Unexpected boundary.');
    end

    E = Erecip + Esurf;

    parts = struct();
    parts.real = zeros(nSites, 3);
    parts.recip = Erecip;
    parts.surf = Esurf;
    parts.self = zeros(nSites, 3);
    parts.total = E;

    parts.Mmu = Mmu;
    parts.Esurf_mu = Esurf_mu;
    parts.surf_coeff = surf_coeff;

    parts.target_mask = targetMask;
    parts.source_mask = sourceMask;
    parts.target_sites = targetSites;
    parts.source_sites = sourceSites;

    parts.nTargets = numel(targetSites);
    parts.nSources = numel(sourceSites);
    parts.nKmesh = solveInfo.nK;

    parts.alpha = alpha;
    parts.boundary = boundary;
    parts.mesh_size = meshSize;
    parts.assignment_order = order;
    parts.deconvolve_assignment = solveInfo.deconvolve_assignment;
    parts.deconvolution_floor = solveInfo.deconvolution_floor;
    parts.derivative_mode = solveInfo.derivative_mode;
    parts.influence_mode = solveInfo.influence_mode;

    parts.use_kcut_mask = solveInfo.use_kcut_mask;
    parts.kcut = solveInfo.kcut;
    parts.prefactor_convention = solveInfo.prefactor_convention;
    parts.nKmesh_total_nonzero = solveInfo.nK_total_nonzero;

    parts.lattice_convention = 'H_rows_G_columns_HG_2piI';
    parts.H = lat.H;
    parts.G = lat.G;
    parts.volume = lat.volume;

    parts.time_mesh = timeMesh;
    parts.solveInfo = solveInfo;

    if verbose
        fprintf('  mesh time   = %.6f s\n', timeMesh);
        fprintf('  Mmu         = [%+.8e %+.8e %+.8e]\n', Mmu);
        fprintf('  ||Erecip||  = %.16e\n', norm(Erecip, 'fro'));
        fprintf('  ||Esurf||   = %.16e\n', norm(Esurf, 'fro'));
        fprintf('  ||E||       = %.16e\n', norm(E, 'fro'));
    end
end

function parts = local_empty_parts(sys, targetMask, sourceMask, alpha, boundary)
    nSites = size(sys.site_pos, 1);

    parts = struct();
    parts.real = zeros(nSites, 3);
    parts.recip = zeros(nSites, 3);
    parts.surf = zeros(nSites, 3);
    parts.self = zeros(nSites, 3);
    parts.total = zeros(nSites, 3);

    parts.Mmu = [0.0, 0.0, 0.0];
    parts.Esurf_mu = [0.0, 0.0, 0.0];
    parts.surf_coeff = 0.0;

    parts.target_mask = targetMask;
    parts.source_mask = sourceMask;
    parts.target_sites = find(targetMask);
    parts.source_sites = find(sourceMask);

    parts.nTargets = nnz(targetMask);
    parts.nSources = nnz(sourceMask);
    parts.nKmesh = 0;

    parts.alpha = alpha;
    parts.boundary = boundary;
    parts.mesh_size = [0 0 0];
    parts.assignment_order = NaN;
    parts.deconvolve_assignment = NaN;
    parts.deconvolution_floor = NaN;
    parts.derivative_mode = 'spectral';
    parts.influence_mode = 'ewald';

    parts.lattice_convention = 'H_rows_G_columns_HG_2piI';
    parts.time_mesh = 0.0;
    parts.solveInfo = struct();
end

function mask = local_get_mask(opts, name, defaultMask)
    if isfield(opts, name) && ~isempty(opts.(name))
        mask = logical(opts.(name)(:));
    else
        mask = logical(defaultMask(:));
    end
end

function x = local_require_struct_field(s, name)
    if ~isfield(s, name) || isempty(s.(name))
        error('p3m:compute_induced_field_dipoles:MissingEwaldField', ...
            'opts.ewald.%s is required.', name);
    end
    x = s.(name);
end

function val = local_get_opt(s, name, defaultVal)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = defaultVal;
    end
end