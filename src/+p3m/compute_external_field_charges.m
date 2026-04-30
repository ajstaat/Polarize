function [E, parts] = compute_external_field_charges(sys, opts)
%COMPUTE_EXTERNAL_FIELD_CHARGES P3M periodic external field from fixed charges.
%
% [E, parts] = p3m.compute_external_field_charges(sys, opts)
%
% Purpose
% -------
% Compute the periodic external electric field from fixed charges using:
%
%   E = E_real,Ewald(+Thole correction) + E_recip,P3M + E_surface
%
% This is the P3M analogue of:
%
%   thole.induced_field_from_charges_periodic
%
% but with the reciprocal charge field evaluated on a mesh.
%
% Project lattice convention
% --------------------------
% Polarize uses direct lattice vectors as ROWS:
%
%   cart = frac * H
%
% geom.get_lattice returns:
%
%   lat.H = H
%   lat.G such that H * G = 2*pi*I
%
% This routine treats geom.get_lattice(sys) as the convention boundary.
%
% Required opts fields
%   opts.ewald.alpha
%   opts.ewald.rcut
%   opts.ewald.kcut       accepted/stored for API symmetry; the P3M mesh
%                         reciprocal solve does not use a spherical kcut
%   opts.mesh_size        [M1 M2 M3]
%   opts.assignment_order B-spline order
%
% Optional opts fields
%   opts.ewald.boundary           'tinfoil' or 'vacuum', default 'tinfoil'
%   opts.target_mask              default true(N,1)
%   opts.source_mask              default abs(site_charge)>0
%   opts.exclude_self             default true
%   opts.use_thole_damping        default true
%   opts.realspace_backend        default 'thole_periodic_real'
%   opts.derivative_mode          default 'spectral'
%   opts.influence_mode           default 'ewald'
%   opts.deconvolve_assignment    default true
%   opts.deconvolution_floor      default 1e-8
%   opts.fd_stencil               default 'central2'
%   opts.alias_range              default 2
%   opts.verbose                  default false
%
% Output
%   E       N x 3 field on all sites, zero on non-target rows
%
%   parts   struct with real/recip/surf pieces and diagnostics

if nargin < 2 || isempty(opts)
    opts = struct();
end

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
rcut = local_require_struct_field(opts.ewald, 'rcut');

% kcut is not used by the mesh reciprocal solve, but keeping it in the
% Ewald parameter struct makes comparison workflows easier and mirrors the
% pure-Ewald external-field API.
kcut = NaN;
if isfield(opts.ewald, 'kcut') && ~isempty(opts.ewald.kcut)
    kcut = opts.ewald.kcut;
end

validateattributes(alpha, {'numeric'}, ...
    {'scalar','real','finite','positive'}, ...
    mfilename, 'opts.ewald.alpha');

validateattributes(rcut, {'numeric'}, ...
    {'scalar','real','finite','positive'}, ...
    mfilename, 'opts.ewald.rcut');

if ~isnan(kcut)
    validateattributes(kcut, {'numeric'}, ...
        {'scalar','real','finite','positive'}, ...
        mfilename, 'opts.ewald.kcut');
end

alpha = double(alpha);
rcut = double(rcut);
kcut = double(kcut);

boundary = 'tinfoil';
if isfield(opts.ewald, 'boundary') && ~isempty(opts.ewald.boundary)
    boundary = lower(char(string(opts.ewald.boundary)));
end

if ~ismember(boundary, {'tinfoil', 'vacuum'})
    error('p3m:compute_external_field_charges:BadBoundary', ...
        'opts.ewald.boundary must be ''tinfoil'' or ''vacuum''.');
end

meshSize = local_validate_mesh_size(opts.mesh_size);
order = local_validate_order(opts.assignment_order);

verbose = local_get_opt(opts, 'verbose', false);
verbose = logical(verbose);

excludeSelf = local_get_opt(opts, 'exclude_self', true);
excludeSelf = logical(excludeSelf);

useTholeDamping = local_get_opt(opts, 'use_thole_damping', true);
useTholeDamping = logical(useTholeDamping);

% Prefer the convention-clean, tested periodic charge-field real cache.
% The old 'matlab' backend is kept only as a fallback for debugging.
realspaceBackend = lower(char(string(local_get_opt( ...
    opts, 'realspace_backend', 'thole_periodic_real'))));

pos = sys.site_pos;
validateattributes(pos, {'numeric'}, ...
    {'2d','ncols',3,'real','finite'}, ...
    mfilename, 'sys.site_pos');

q = sys.site_charge(:);
nSites = size(pos, 1);

if numel(q) ~= nSites
    error('p3m:compute_external_field_charges:BadChargeSize', ...
        'sys.site_charge must have length N.');
end

targetMask = local_get_mask(opts, 'target_mask', true(nSites, 1));
sourceMask = local_get_mask(opts, 'source_mask', abs(q) > 0);

if numel(targetMask) ~= nSites || numel(sourceMask) ~= nSites
    error('p3m:compute_external_field_charges:BadMaskSize', ...
        'target_mask and source_mask must have length N.');
end

targetMask = logical(targetMask(:));
sourceMask = logical(sourceMask(:));

targetSites = find(targetMask);
sourceSites = find(sourceMask);

lat = geom.get_lattice(sys);
H = lat.H;
V = lat.volume;

if V <= 1e-14
    error('p3m:compute_external_field_charges:BadVolume', ...
        'Cell volume is nonpositive.');
end

qSource = q(sourceSites);
qtot = sum(qSource);

if abs(qtot) > 1e-10
    error('p3m:compute_external_field_charges:NonNeutralSources', ...
        ['Periodic P3M field from fixed charges requires net-neutral ', ...
         'selected sources. Selected total charge = %+0.16e'], qtot);
end

if verbose
    fprintf('p3m.compute_external_field_charges:\n');
    fprintf('  lattice convention       = project_row_H_column_G_HG_2piI\n');
    fprintf('  mesh_size                = [%d %d %d]\n', meshSize);
    fprintf('  assignment_order         = %d\n', order);
    fprintf('  alpha                    = %.8g\n', alpha);
    fprintf('  rcut                     = %.8g bohr\n', rcut);
    if ~isnan(kcut)
        fprintf('  kcut                     = %.8g bohr^-1 metadata only\n', kcut);
    end
    fprintf('  boundary                 = %s\n', boundary);
    fprintf('  realspace_backend        = %s\n', realspaceBackend);
    fprintf('  use_thole_damping        = %d\n', useTholeDamping);
    fprintf('  exclude_self             = %d\n', excludeSelf);
    fprintf('  nTargets                 = %d\n', numel(targetSites));
    fprintf('  nSources                 = %d\n', numel(sourceSites));
end

%% ------------------------------------------------------------------------
% Real-space Ewald / short-range correction piece
% -------------------------------------------------------------------------

tReal = tic;

switch realspaceBackend
    case 'thole_periodic_real'
        fieldReal = struct();
        fieldReal.mode = 'periodic';
        fieldReal.real_only = true;
        fieldReal.exclude_self = excludeSelf;
        fieldReal.use_thole_damping = useTholeDamping;
        fieldReal.target_mask = targetMask;
        fieldReal.source_mask = sourceMask;
        fieldReal.verbose = false;

        % These are accepted by the periodic field driver even though the
        % real-only path should not need reciprocal storage.
        fieldReal.kspace_mode = 'full';
        fieldReal.k_block_size = local_get_opt(opts, 'k_block_size', 2048);
        fieldReal.kspace_memory_limit_gb = local_get_opt(opts, ...
            'kspace_memory_limit_gb', 8);

        fieldReal.ewald = opts.ewald;
        fieldReal.ewald.alpha = alpha;
        fieldReal.ewald.rcut = rcut;
        fieldReal.ewald.boundary = boundary;

        [~, realParts] = thole.induced_field_from_charges_periodic(sys, fieldReal);

        Ereal = realParts.real;

        realInfo = struct();
        realInfo.backend = 'thole_periodic_real';
        realInfo.nTerms = local_get_opt(realParts, 'nRealEntries', NaN);
        realInfo.alpha = alpha;
        realInfo.rcut = rcut;
        realInfo.time_real = local_get_opt(realParts, 'time_real', NaN);
        realInfo.real_only = local_get_opt(realParts, 'real_only', true);

        if isfield(realParts, 'realCache')
            realInfo.realCache = realParts.realCache;
        end

    case 'matlab'
        % Legacy/debug path only. Keep available if old local P3M files
        % include p3m.realspace_charge_field_ewald.
        realOpts = struct();
        realOpts.alpha = alpha;
        realOpts.rcut = rcut;
        realOpts.target_mask = targetMask;
        realOpts.source_mask = sourceMask;
        realOpts.exclude_self = excludeSelf;
        realOpts.use_thole_damping = useTholeDamping;

        if exist('p3m.realspace_charge_field_ewald', 'file') ~= 2
            error('p3m:compute_external_field_charges:MissingLegacyRealspace', ...
                ['realspace_backend=''matlab'' requested, but ', ...
                 'p3m.realspace_charge_field_ewald is not available.']);
        end

        [Ereal, realInfo] = p3m.realspace_charge_field_ewald(sys, realOpts);

    otherwise
        error('p3m:compute_external_field_charges:BadRealspaceBackend', ...
            'Unknown realspace_backend: %s', realspaceBackend);
end

timeReal = toc(tReal);

%% ------------------------------------------------------------------------
% Mesh reciprocal P3M piece
% -------------------------------------------------------------------------

tMesh = tic;

fracSource = pos(sourceSites, :) / H;
fracTarget = pos(targetSites, :) / H;

rho = p3m.assign_charges_bsplines(fracSource, qSource, meshSize, order);

solveOpts = struct();
solveOpts.alpha = alpha;
solveOpts.assignment_order = order;
solveOpts.deconvolve_assignment = local_get_opt(opts, ...
    'deconvolve_assignment', true);
solveOpts.deconvolution_floor = local_get_opt(opts, ...
    'deconvolution_floor', 1e-8);
solveOpts.derivative_mode = local_get_opt(opts, ...
    'derivative_mode', 'spectral');
solveOpts.fd_stencil = local_get_opt(opts, ...
    'fd_stencil', 'central2');
solveOpts.influence_mode = local_get_opt(opts, ...
    'influence_mode', 'ewald');
solveOpts.alias_range = local_get_opt(opts, ...
    'alias_range', 2);

[ExGrid, EyGrid, EzGrid, solveInfo] = ...
    p3m.solve_charge_field_spectral(rho, lat, solveOpts);

ErecipTarget = p3m.interpolate_field_bsplines( ...
    fracTarget, ExGrid, EyGrid, EzGrid, order);

Erecip = zeros(nSites, 3);
Erecip(targetSites, :) = ErecipTarget;

timeMesh = toc(tMesh);

%% ------------------------------------------------------------------------
% Surface term
% -------------------------------------------------------------------------

Esurf = zeros(nSites, 3);

% Origin-dependent for a neutral source distribution, as usual; this is the
% same convention used in the pure Ewald fixed-charge field path.
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
        error('p3m:compute_external_field_charges:BadBoundaryInternal', ...
            'Unexpected boundary "%s".', boundary);
end

%% ------------------------------------------------------------------------
% Pack result
% -------------------------------------------------------------------------

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

parts.alpha = alpha;
parts.rcut = rcut;
parts.kcut = kcut;
parts.boundary = boundary;

parts.realspace_backend = realspaceBackend;
parts.use_thole_damping = useTholeDamping;
parts.exclude_self = excludeSelf;

parts.derivative_mode = solveInfo.derivative_mode;
parts.fd_stencil = solveInfo.fd_stencil;
parts.influence_mode = solveInfo.influence_mode;
parts.alias_range = solveInfo.alias_range;

parts.lattice = lat;
parts.H = lat.H;
parts.G = lat.G;
parts.volume = lat.volume;
parts.lattice_convention = 'project_row_H_column_G_HG_2piI';

parts.realInfo = realInfo;
parts.solveInfo = solveInfo;

parts.time_real = timeReal;
parts.time_mesh = timeMesh;
parts.time_total = timeReal + timeMesh;

parts.rho_total = sum(rho(:));
parts.rho = rho;

parts.grid_field = struct();
parts.grid_field.Ex = ExGrid;
parts.grid_field.Ey = EyGrid;
parts.grid_field.Ez = EzGrid;

parts.norm_Ereal = norm(Ereal, 'fro');
parts.norm_Erecip = norm(Erecip, 'fro');
parts.norm_Esurf = norm(Esurf, 'fro');
parts.norm_Etotal = norm(E, 'fro');

if isfield(solveInfo, 'nK')
    parts.nK = solveInfo.nK;
else
    parts.nK = NaN;
end

if verbose
    fprintf('  real time                = %.6f s\n', timeReal);
    fprintf('  mesh time                = %.6f s\n', timeMesh);
    fprintf('  total time               = %.6f s\n', parts.time_total);
    fprintf('  sum q source             = %+0.16e\n', qtot);
    fprintf('  sum rho mesh             = %+0.16e\n', parts.rho_total);
    fprintf('  ||Ereal||_F              = %.16e\n', parts.norm_Ereal);
    fprintf('  ||Erecip||_F             = %.16e\n', parts.norm_Erecip);
    fprintf('  ||Esurf||_F              = %.16e\n', parts.norm_Esurf);
    fprintf('  ||E||_F                  = %.16e\n', parts.norm_Etotal);
end
end

% =========================================================================
% Local helpers
% =========================================================================

function meshSize = local_validate_mesh_size(meshSize)
validateattributes(meshSize, {'numeric'}, ...
    {'vector','numel',3,'integer','positive','finite'}, ...
    mfilename, 'opts.mesh_size');

meshSize = double(meshSize(:).');
end

function order = local_validate_order(order)
validateattributes(order, {'numeric'}, ...
    {'scalar','integer','positive','finite'}, ...
    mfilename, 'opts.assignment_order');

order = double(order);
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