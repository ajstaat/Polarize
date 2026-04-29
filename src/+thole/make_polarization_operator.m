function op = make_polarization_operator(sys, problem, varargin)
%MAKE_POLARIZATION_OPERATOR Build a polarization operator object.
%
% op = thole.make_polarization_operator(sys, problem, Name, Value)
%
% Public options
%   'Mode'
%       'nonperiodic', 'periodic_ewald', 'periodic_p3m'
%       default: inferred from sys.is_periodic
%
%   'Solver'
%       'direct', 'jacobi', 'gmres', 'sor'
%       default: 'direct'
%
%   'Backend'
%       'auto', 'dense', 'matrix_free'
%       default: 'auto'
%
%   'UseThole'
%       logical, default true
%
%   'Softening'
%       scalar nonnegative softening, default 0
%       Periodic Ewald operators currently require Softening = 0.
%
%   'Rcut'
%       real-space cutoff in bohr.
%       Required finite for matrix-free nonperiodic and all periodic Ewald
%       operators.
%
%   'UseMex'
%       logical, default true
%
%   'Profile'
%       logical, default false
%
%   'Verbose'
%       logical, default false
%
% Periodic Ewald options
%   These may be supplied either as top-level Name/Value pairs or through
%   'Ewald', struct(...).
%
%   'Alpha' or Ewald.alpha
%       positive Ewald screening parameter
%
%   'Kcut' or Ewald.kcut
%       positive reciprocal-space cutoff
%
%   'Boundary' or Ewald.boundary
%       'tinfoil' or 'vacuum', default 'tinfoil'
%
%   'KspaceMode'
%       'auto', 'full', 'blocked', default 'auto'
%
%   'KBlockSize'
%       positive integer, default 2048
%
%   'KspaceMemoryLimitGB'
%       positive scalar, default 8
%
%   'UseMexKspace'
%       logical, default true
%
% Operator fields
%   op.kind
%       'dense_matrix' or 'matrix_free'
%
%   op.backend
%       Concrete implementation, e.g.
%         'nonperiodic_paircache_dense'
%         'nonperiodic_paircache_apply'
%         'nonperiodic_rowcache_apply'
%         'periodic_paircache_dense'
%         'periodic_paircache_apply'
%         'periodic_rowcache_apply'

p = inputParser;

addRequired(p, 'sys', @isstruct);
addRequired(p, 'problem', @isstruct);

addParameter(p, 'Mode', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Solver', 'direct', @(x) ischar(x) || isstring(x));
addParameter(p, 'Backend', 'auto', @(x) ischar(x) || isstring(x));

addParameter(p, 'UseThole', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Softening', 0.0, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);

addParameter(p, 'Rcut', Inf, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(p, 'UseMex', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Profile', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));

% Periodic Ewald options. These are optional at parser level because
% nonperiodic operators do not need them. Periodic validation happens after
% mode resolution.
addParameter(p, 'Ewald', struct(), @(x) isempty(x) || isstruct(x));

addParameter(p, 'Alpha', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));

addParameter(p, 'Kcut', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));

addParameter(p, 'Boundary', 'tinfoil', @(x) ischar(x) || isstring(x));

addParameter(p, 'KspaceMode', 'auto', @(x) ischar(x) || isstring(x));

addParameter(p, 'KBlockSize', 2048, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x == round(x));

addParameter(p, 'KspaceMemoryLimitGB', 8, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);

addParameter(p, 'UseMexKspace', true, @(x) islogical(x) && isscalar(x));

parse(p, sys, problem, varargin{:});
opt = p.Results;

mode = local_normalize_mode(opt.Mode, sys);
solver = local_normalize_solver(opt.Solver);
backendRequest = local_normalize_backend_request(opt.Backend);

validate_problem(problem);

if strcmp(mode, 'periodic_p3m')
    error('thole:make_polarization_operator:PeriodicP3MNotImplemented', ...
        ['Mode="periodic_p3m" is recognized, but P3M operators are not ', ...
         'implemented in this refactor stage. Use Mode="periodic_ewald".']);
end

if strcmp(mode, 'periodic_ewald')
    opt = local_normalize_periodic_ewald_options(opt);
end

backendResolved = local_resolve_backend(mode, solver, backendRequest, opt.Rcut);

switch backendResolved
    case 'dense'
        op = build_nonperiodic_dense_operator(sys, problem, opt);

    case 'matrix_free_paircache'
        op = build_nonperiodic_paircache_operator(sys, problem, opt);

    case 'matrix_free_rowcache'
        op = build_nonperiodic_rowcache_operator(sys, problem, opt);

    case 'periodic_dense'
        op = build_periodic_dense_operator(sys, problem, opt);

    case 'periodic_paircache'
        op = build_periodic_paircache_operator(sys, problem, opt);

    case 'periodic_rowcache'
        op = build_periodic_rowcache_operator(sys, problem, opt);

    otherwise
        error('thole:make_polarization_operator:InternalBadBackendResolution', ...
            'Internal error: unresolved backend "%s".', backendResolved);
end

op.request = struct();
op.request.mode = mode;
op.request.solver = solver;
op.request.backend = backendRequest;
op.request.backend_resolved = backendResolved;
end

% =========================================================================
% Resolution / validation
% =========================================================================

function mode = local_normalize_mode(modeIn, sys)
mode = lower(strtrim(char(string(modeIn))));

if isempty(mode)
    if isfield(sys, 'is_periodic') && ~isempty(sys.is_periodic) && logical(sys.is_periodic)
        mode = 'periodic_ewald';
    else
        mode = 'nonperiodic';
    end
end

switch mode
    case {'nonperiodic', 'finite', 'cluster'}
        mode = 'nonperiodic';

    case {'periodic', 'periodic_ewald', 'ewald'}
        mode = 'periodic_ewald';

    case {'periodic_p3m', 'p3m'}
        mode = 'periodic_p3m';

    otherwise
        error('thole:make_polarization_operator:UnsupportedOperatorMode', ...
            'Unsupported polarization-operator mode "%s".', mode);
end
end

function solver = local_normalize_solver(solverIn)
solver = lower(strtrim(char(string(solverIn))));

switch solver
    case {'direct', 'dense_direct'}
        solver = 'direct';

    case {'jacobi', 'iterative', 'fixed_point'}
        solver = 'jacobi';

    case {'gmres', 'krylov'}
        solver = 'gmres';

    case {'sor', 'gauss_seidel'}
        solver = 'sor';

    otherwise
        error('thole:make_polarization_operator:UnsupportedSolver', ...
            ['Unsupported solver "%s". Supported solver intents are: ', ...
             '"direct", "jacobi", "gmres", "sor".'], solver);
end
end

function backend = local_normalize_backend_request(backendIn)
backend = lower(strtrim(char(string(backendIn))));

switch backend
    case {'auto', ''}
        backend = 'auto';

    case {'dense', 'dense_matrix'}
        backend = 'dense';

    case {'matrix_free', 'matrixfree'}
        backend = 'matrix_free';

    case {'pair_cache', 'paircache'}
        error('thole:make_polarization_operator:PairCacheBackendRenamed', ...
            ['Backend="pair_cache" has been replaced by Backend="matrix_free". ', ...
             'The pair cache is an internal implementation detail.']);

    case {'row_cache', 'rowcache'}
        error('thole:make_polarization_operator:RowCacheBackendIsInternal', ...
            ['Backend="row_cache" is an internal implementation detail. ', ...
             'Use Backend="matrix_free" with Solver="sor".']);

    otherwise
        error('thole:make_polarization_operator:UnsupportedOperatorBackend', ...
            ['Unsupported Backend="%s". Supported public backend values are: ', ...
             '"auto", "dense", "matrix_free".'], backend);
end
end

function backend = local_resolve_backend(mode, solver, backendRequest, rcut)
switch mode
    case 'nonperiodic'
        backend = local_resolve_nonperiodic_backend(solver, backendRequest, rcut);

    case 'periodic_ewald'
        backend = local_resolve_periodic_ewald_backend(solver, backendRequest, rcut);

    otherwise
        error('thole:make_polarization_operator:UnsupportedOperatorMode', ...
            'Unsupported operator mode "%s".', mode);
end
end

function backend = local_resolve_nonperiodic_backend(solver, backendRequest, rcut)
if strcmp(backendRequest, 'dense')
    if strcmp(solver, 'sor')
        error('thole:make_polarization_operator:SorRequiresMatrixFreeRowUpdate', ...
            ['SOR requires a matrix-free operator with row-update capability. ', ...
             'Backend="dense" is not compatible with SOR.']);
    end

    backend = 'dense';
    return;
end

if strcmp(backendRequest, 'matrix_free')
    if strcmp(solver, 'direct')
        error('thole:make_polarization_operator:DirectRequiresDenseOperator', ...
            ['Direct SCF requires a dense operator with op.Tpol. ', ...
             'Use Backend="dense" or Backend="auto" with Solver="direct".']);
    end

    if ~isfinite(rcut)
        error('thole:make_polarization_operator:MatrixFreeRequiresFiniteCutoff', ...
            ['Nonperiodic matrix-free operator construction currently requires ', ...
             'a finite Rcut. Use Backend="dense" for full all-pairs calculations, ', ...
             'or provide a finite Rcut.']);
    end

    if strcmp(solver, 'sor')
        backend = 'matrix_free_rowcache';
    else
        backend = 'matrix_free_paircache';
    end

    return;
end

% Backend = auto
switch solver
    case 'direct'
        backend = 'dense';

    case {'jacobi', 'gmres'}
        if ~isfinite(rcut)
            error('thole:make_polarization_operator:AutoMatrixFreeRequiresFiniteCutoff', ...
                ['Backend="auto" with Solver="%s" resolves to a matrix-free ', ...
                 'operator, which currently requires a finite Rcut for ', ...
                 'nonperiodic systems. Use Backend="dense" for full all-pairs ', ...
                 'calculations, or provide a finite Rcut.'], solver);
        end

        backend = 'matrix_free_paircache';

    case 'sor'
        if ~isfinite(rcut)
            error('thole:make_polarization_operator:SorRequiresFiniteCutoff', ...
                ['Backend="auto" with Solver="sor" resolves to a matrix-free ', ...
                 'row-cache operator, which currently requires a finite Rcut ', ...
                 'for nonperiodic systems. Provide a finite Rcut.']);
        end

        backend = 'matrix_free_rowcache';

    otherwise
        error('thole:make_polarization_operator:UnsupportedSolver', ...
            'Unsupported solver "%s".', solver);
end
end

function backend = local_resolve_periodic_ewald_backend(solver, backendRequest, rcut)
if ~isfinite(rcut)
    error('thole:make_polarization_operator:PeriodicEwaldRequiresFiniteRcut', ...
        ['Periodic Ewald operators require a finite real-space Rcut. ', ...
         'Provide Rcut or Ewald.rcut.']);
end

if strcmp(backendRequest, 'dense')
    if strcmp(solver, 'sor')
        error('thole:make_polarization_operator:SorRequiresMatrixFreeRowUpdate', ...
            ['SOR requires a matrix-free periodic row-cache operator with ', ...
             'row-update capability. Backend="dense" is not compatible with SOR.']);
    end

    backend = 'periodic_dense';
    return;
end

if strcmp(backendRequest, 'matrix_free')
    if strcmp(solver, 'direct')
        error('thole:make_polarization_operator:DirectRequiresDenseOperator', ...
            ['Direct SCF requires a dense operator with op.Tpol. ', ...
             'Use Backend="dense" or Backend="auto" with Solver="direct".']);
    end

    if strcmp(solver, 'sor')
        backend = 'periodic_rowcache';
    else
        backend = 'periodic_paircache';
    end

    return;
end

% Backend = auto
switch solver
    case 'direct'
        backend = 'periodic_dense';

    case {'jacobi', 'gmres'}
        backend = 'periodic_paircache';

    case 'sor'
        backend = 'periodic_rowcache';

    otherwise
        error('thole:make_polarization_operator:UnsupportedSolver', ...
            'Unsupported solver "%s".', solver);
end
end

function opt = local_normalize_periodic_ewald_options(opt)
if isfield(opt, 'Ewald') && ~isempty(opt.Ewald)
    ew = opt.Ewald;
else
    ew = struct();
end

alpha = local_get_first_field(opt, ew, {'Alpha', 'alpha'}, []);
kcut = local_get_first_field(opt, ew, {'Kcut', 'kcut'}, []);
rcut = local_get_first_field(opt, ew, {'Rcut', 'rcut'}, opt.Rcut);

if isempty(alpha)
    error('thole:make_polarization_operator:MissingPeriodicAlpha', ...
        'Mode="periodic_ewald" requires Alpha or Ewald.alpha.');
end

if isempty(kcut)
    error('thole:make_polarization_operator:MissingPeriodicKcut', ...
        'Mode="periodic_ewald" requires Kcut or Ewald.kcut.');
end

if isempty(rcut) || ~isfinite(rcut)
    error('thole:make_polarization_operator:MissingPeriodicRcut', ...
        'Mode="periodic_ewald" requires finite Rcut or Ewald.rcut.');
end

validateattributes(alpha, {'numeric'}, ...
    {'scalar','real','finite','positive'}, ...
    mfilename, 'Alpha');

validateattributes(kcut, {'numeric'}, ...
    {'scalar','real','finite','positive'}, ...
    mfilename, 'Kcut');

validateattributes(rcut, {'numeric'}, ...
    {'scalar','real','finite','positive'}, ...
    mfilename, 'Rcut');

boundary = local_get_first_field(opt, ew, {'Boundary', 'boundary'}, opt.Boundary);
boundary = lower(char(string(boundary)));

if ~ismember(boundary, {'tinfoil','vacuum'})
    error('thole:make_polarization_operator:BadPeriodicBoundary', ...
        'Periodic Ewald Boundary must be ''tinfoil'' or ''vacuum''.');
end

kspaceMode = lower(char(string(opt.KspaceMode)));
if strcmp(kspaceMode, 'chunked')
    kspaceMode = 'blocked';
end

if ~ismember(kspaceMode, {'auto','full','blocked'})
    error('thole:make_polarization_operator:BadKspaceMode', ...
        'KspaceMode must be ''auto'', ''full'', or ''blocked''.');
end

if opt.Softening ~= 0
    error('thole:make_polarization_operator:PeriodicSofteningUnsupported', ...
        'Mode="periodic_ewald" currently requires Softening = 0.');
end

opt.Alpha = double(alpha);
opt.Kcut = double(kcut);
opt.Rcut = double(rcut);
opt.Boundary = boundary;
opt.KspaceMode = kspaceMode;
opt.KBlockSize = double(opt.KBlockSize);
opt.KspaceMemoryLimitGB = double(opt.KspaceMemoryLimitGB);

opt.Ewald = struct();
opt.Ewald.alpha = opt.Alpha;
opt.Ewald.kcut = opt.Kcut;
opt.Ewald.rcut = opt.Rcut;
opt.Ewald.boundary = opt.Boundary;
end

function value = local_get_first_field(opt, ew, names, defaultValue)
value = defaultValue;

for k = 1:numel(names)
    name = names{k};

    if isfield(opt, name) && ~isempty(opt.(name))
        value = opt.(name);
        return;
    end
end

for k = 1:numel(names)
    name = names{k};

    if isfield(ew, name) && ~isempty(ew.(name))
        value = ew.(name);
        return;
    end
end
end

function validate_problem(problem)
required = {'activeSites', 'nPolSites', 'Eext_pol_vec', 'alpha_pol_vec'};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(problem, name) || isempty(problem.(name))
        error('thole:make_polarization_operator:MissingProblemField', ...
            'problem.%s is required and missing/empty.', name);
    end
end

nVec = 3 * problem.nPolSites;

if numel(problem.Eext_pol_vec) ~= nVec
    error('thole:make_polarization_operator:BadProblemSize', ...
        'problem.Eext_pol_vec must have length 3*problem.nPolSites.');
end

if numel(problem.alpha_pol_vec) ~= nVec
    error('thole:make_polarization_operator:BadProblemSize', ...
        'problem.alpha_pol_vec must have length 3*problem.nPolSites.');
end

if numel(problem.activeSites) ~= problem.nPolSites
    error('thole:make_polarization_operator:BadProblemSize', ...
        'numel(problem.activeSites) must equal problem.nPolSites.');
end
end