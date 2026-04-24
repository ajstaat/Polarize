function [E, parts] = compute_external_field_charges(sys, opts)
%COMPUTE_EXTERNAL_FIELD_CHARGES P3M periodic external field from charges.
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
% This routine uses geom.get_lattice(sys) and passes the resulting lattice
% struct into P3M reciprocal solvers.
%
% Computes:
%
%   E = E_real + E_recip,P3M + E_surface

    io.assert_atomic_units(sys);

    if ~isfield(opts, 'ewald') || isempty(opts.ewald)
        error('p3m:compute_external_field_charges:MissingEwald', ...
            'opts.ewald is required.');
    end

    if ~isfield(opts, 'mesh_size') || isempty(opts.mesh_size)
        error('p3m:compute_external_field_charges:MissingMeshSize', ...
            'opts.mesh_size is required.');
    end

    if ~isfield(opts, 'assignment_order') || isempty(opts.assignment_order)
        error('p3m:compute_external_field_charges:MissingAssignmentOrder', ...
            'opts.assignment_order is required.');
    end

    alpha = local_require_struct_field(opts.ewald, 'alpha');
    rcut  = local_require_struct_field(opts.ewald, 'rcut');

    boundary = 'tinfoil';
    if isfield(opts.ewald, 'boundary') && ~isempty(opts.ewald.boundary)
        boundary = lower(char(string(opts.ewald.boundary)));
    end
    if ~ismember(boundary, {'tinfoil', 'vacuum'})
        error('p3m:compute_external_field_charges:BadBoundary', ...
            'boundary must be ''tinfoil'' or ''vacuum''.');
    end

    meshSize = double(opts.mesh_size(:).');
    if numel(meshSize) ~= 3 || any(meshSize < 1) || any(meshSize ~= round(meshSize))
        error('p3m:compute_external_field_charges:BadMeshSize', ...
            'opts.mesh_size must be positive integer [M1 M2 M3].');
    end

    order = opts.assignment_order;
    validateattributes(order, {'numeric'}, ...
        {'scalar','real','finite','positive','integer'}, ...
        mfilename, 'opts.assignment_order');
    order = double(order);

    verbose = local_get_opt(opts, 'verbose', false);
    realspaceBackend = lower(char(string(local_get_opt(opts, 'realspace_backend', 'matlab'))));
    useTholeDamping = local_get_opt(opts, 'use_thole_damping', false);

    pos = sys.site_pos;
    q = sys.site_charge(:);
    nSites = size(pos, 1);

    if numel(q) ~= nSites
        error('p3m:compute_external_field_charges:BadChargeSize', ...
            'sys.site_charge must have length N.');
    end

    targetMask = local_get_mask(opts, 'target_mask', true(nSites,1));
    sourceMask = local_get_mask(opts, 'source_mask', abs(q) > 0);

    if numel(targetMask) ~= nSites || numel(sourceMask) ~= nSites
        error('p3m:compute_external_field_charges:BadMask', ...
            'target_mask and source_mask must have length N.');
    end

    lat = geom.get_lattice(sys);
    H = lat.H;
    V = lat.volume;

    if V <= 1e-14
        error('p3m:compute_external_field_charges:BadVolume', ...
            'Cell volume is nonpositive.');
    end

    if verbose
        fprintf('p3m.compute_external_field_charges:\n');
        fprintf('  lattice convention = H_rows_G_columns_HG_2piI\n');
        fprintf('  mesh_size         = [%d %d %d]\n', meshSize);
        fprintf('  order             = %d\n', order);
        fprintf('  alpha             = %.8g\n', alpha);
        fprintf('  rcut              = %.8g\n', rcut);
        fprintf('  boundary          = %s\n', boundary);
        fprintf('  realspace_backend = %s\n', realspaceBackend);
        fprintf('  use_thole_damping = %d\n', useTholeDamping);
        fprintf('  nTargets          = %d\n', nnz(targetMask));
        fprintf('  nSources          = %d\n', nnz(sourceMask));
    end

    %% --------------------------------------------------------------------
    % Real-space Ewald part

    tReal = tic;

    switch realspaceBackend
        case 'matlab'
            realOpts = struct();
            realOpts.alpha = alpha;
            realOpts.rcut = rcut;
            realOpts.target_mask = targetMask;
            realOpts.source_mask = sourceMask;
            realOpts.exclude_self = local_get_opt(opts, 'exclude_self', true);

            [Ereal, realInfo] = p3m.realspace_charge_field_ewald(sys, realOpts);

        case 'thole_periodic_real'
            fieldReal = struct();
            fieldReal.mode = 'periodic';
            fieldReal.real_only = true;
            fieldReal.exclude_self = local_get_opt(opts, 'exclude_self', true);
            fieldReal.use_thole_damping = useTholeDamping;
            fieldReal.target_mask = targetMask;
            fieldReal.source_mask = sourceMask;
            fieldReal.verbose = false;

            fieldReal.kspace_mode = 'full';
            fieldReal.k_block_size = local_get_opt(opts, 'k_block_size', 2048);
            fieldReal.kspace_memory_limit_gb = local_get_opt(opts, 'kspace_memory_limit_gb', 8);
            fieldReal.ewald = opts.ewald;

            [~, realParts] = thole.induced_field_from_charges_periodic(sys, fieldReal);

            Ereal = realParts.real;

            realInfo = struct();
            realInfo.backend = 'thole_periodic_real';
            realInfo.nTerms = local_get_opt(realParts, 'nRealEntries', NaN);
            realInfo.alpha = alpha;
            realInfo.rcut = rcut;
            realInfo.time_real = local_get_opt(realParts, 'time_real', NaN);
            realInfo.real_only = local_get_opt(realParts, 'real_only', true);

        otherwise
            error('p3m:compute_external_field_charges:BadRealspaceBackend', ...
                'Unknown realspace_backend: %s', realspaceBackend);
    end

    timeReal = toc(tReal);

    %% --------------------------------------------------------------------
    % Mesh reciprocal part

    sourceSites = find(sourceMask);
    targetSites = find(targetMask);

    fracSource = pos(sourceSites, :) / H;
    fracTarget = pos(targetSites, :) / H;

    qSource = q(sourceSites);

    qtot = sum(qSource);
    if abs(qtot) > 1e-10
        error('p3m:compute_external_field_charges:NonNeutralSources', ...
            ['Periodic P3M field from fixed charges requires net-neutral selected sources. ' ...
             'Selected total charge = %+0.16e'], qtot);
    end

    tMesh = tic;

    rho = p3m.assign_charges_bsplines(fracSource, qSource, meshSize, order);

    solveOpts = struct();
    solveOpts.alpha = alpha;
    solveOpts.assignment_order = order;
    solveOpts.deconvolve_assignment = local_get_opt(opts, 'deconvolve_assignment', true);
    solveOpts.deconvolution_floor = local_get_opt(opts, 'deconvolution_floor', 1e-8);
    solveOpts.derivative_mode = local_get_opt(opts, 'derivative_mode', 'spectral');
    solveOpts.fd_stencil = local_get_opt(opts, 'fd_stencil', 'central2');
    solveOpts.influence_mode = local_get_opt(opts, 'influence_mode', 'ewald');
    solveOpts.alias_range = local_get_opt(opts, 'alias_range', 2);

    [ExGrid, EyGrid, EzGrid, solveInfo] = ...
        p3m.solve_charge_field_spectral(rho, lat, solveOpts);

    ErecipTarget = p3m.interpolate_field_bsplines( ...
        ExGrid, EyGrid, EzGrid, fracTarget, meshSize, order);

    Erecip = zeros(nSites, 3);
    Erecip(targetSites, :) = ErecipTarget;

    timeMesh = toc(tMesh);

    %% --------------------------------------------------------------------
    % Surface term

    Esurf = zeros(nSites, 3);

    Mq = sum(qSource .* pos(sourceSites, :), 1);

    switch boundary
        case 'tinfoil'
            Esurf_q = [0.0, 0.0, 0.0];
            surf_coeff = 0.0;

        case 'vacuum'
            surf_coeff = 4*pi/(3*V);
            Esurf_q = -surf_coeff * Mq;
            Esurf(targetSites, :) = repmat(Esurf_q, numel(targetSites), 1);

        otherwise
            error('p3m:compute_external_field_charges:BadBoundary', ...
                'boundary must be ''tinfoil'' or ''vacuum''.');
    end

    E = Ereal + Erecip + Esurf;

    parts = struct();
    parts.real = Ereal;
    parts.recip = Erecip;
    parts.surf = Esurf;
    parts.total = E;

    parts.Mq = Mq;
    parts.Esurf_q = Esurf_q;
    parts.surf_coeff = surf_coeff;
    parts.qtot = qtot;
    parts.target_mask = targetMask;
    parts.source_mask = sourceMask;
    parts.target_sites = targetSites;
    parts.source_sites = sourceSites;

    parts.mesh_size = meshSize;
    parts.assignment_order = order;
    parts.derivative_mode = solveInfo.derivative_mode;
    parts.fd_stencil = solveInfo.fd_stencil;
    parts.influence_mode = solveInfo.influence_mode;
    parts.alias_range = solveInfo.alias_range;
    parts.alpha = alpha;
    parts.rcut = rcut;
    parts.boundary = boundary;
    parts.realspace_backend = realspaceBackend;
    parts.use_thole_damping = useTholeDamping;

    parts.lattice_convention = 'H_rows_G_columns_HG_2piI';
    parts.H = lat.H;
    parts.G = lat.G;
    parts.volume = lat.volume;

    parts.realInfo = realInfo;
    parts.solveInfo = solveInfo;

    parts.time_real = timeReal;
    parts.time_mesh = timeMesh;
    parts.time_total = timeReal + timeMesh;

    parts.rho_total = sum(rho(:));

    if verbose
        fprintf('  real time  = %.6f s\n', timeReal);
        fprintf('  mesh time  = %.6f s\n', timeMesh);
        fprintf('  total time = %.6f s\n', parts.time_total);
        fprintf('  sum q source = %+0.16e\n', qtot);
        fprintf('  sum rho mesh = %+0.16e\n', parts.rho_total);
        fprintf('  ||Ereal||_F  = %.16e\n', norm(Ereal, 'fro'));
        fprintf('  ||Erecip||_F = %.16e\n', norm(Erecip, 'fro'));
        fprintf('  ||Esurf||_F  = %.16e\n', norm(Esurf, 'fro'));
        fprintf('  ||E||_F      = %.16e\n', norm(E, 'fro'));
    end
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
        error('p3m:compute_external_field_charges:MissingEwaldField', ...
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