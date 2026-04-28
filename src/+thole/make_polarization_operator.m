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
%   'UseThole'   logical, default true
%   'Softening'  scalar, default 0
%   'Rcut'       scalar cutoff in bohr, default Inf
%   'UseMex'     logical, default true
%   'Profile'    logical, default false
%   'Verbose'    logical, default false
%
% Operator fields
%   op.kind
%       Broad representation class used by solvers:
%           'dense_matrix'
%           'matrix_free'
%
%   op.backend
%       Concrete implementation:
%           'nonperiodic_paircache_dense'
%           'nonperiodic_allpairs_dense'
%           'nonperiodic_paircache_apply'
%
% Current support
%   nonperiodic + direct + auto/dense
%   nonperiodic + jacobi/gmres + auto/matrix_free with finite Rcut
%   nonperiodic + jacobi/gmres + dense
%
% SOR row-cache support will be added as another private backend builder.

p = inputParser;
addRequired(p, 'sys', @isstruct);
addRequired(p, 'problem', @isstruct);

addParameter(p, 'Mode', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Solver', 'direct', @(x) ischar(x) || isstring(x));
addParameter(p, 'Backend', 'auto', @(x) ischar(x) || isstring(x));

addParameter(p, 'UseThole', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Softening', 0.0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'Rcut', Inf, @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(p, 'UseMex', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Profile', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));

parse(p, sys, problem, varargin{:});
opt = p.Results;

mode = local_normalize_mode(opt.Mode, sys);
solver = local_normalize_solver(opt.Solver);
backendRequest = local_normalize_backend_request(opt.Backend);

validate_problem(problem);

if ~strcmp(mode, 'nonperiodic')
    error('thole:make_polarization_operator:UnsupportedOperatorMode', ...
        ['Unsupported polarization operator mode "%s" in this refactor stage. ', ...
         'Currently supported: mode="nonperiodic".'], mode);
end

backendResolved = local_resolve_backend(mode, solver, backendRequest, opt.Rcut);

switch backendResolved
    case 'dense'
        op = build_nonperiodic_dense_operator(sys, problem, opt);

    case 'matrix_free_paircache'
        op = build_nonperiodic_paircache_operator(sys, problem, opt);

    case 'matrix_free_rowcache'
        op = build_nonperiodic_rowcache_operator(sys, problem, opt);

    otherwise
        error('thole:make_polarization_operator:InternalBadBackendResolution', ...
            'Internal error: unresolved backend "%s".', backendResolved);
end

op.request = struct();
op.request.mode = mode;
op.request.solver = solver;
op.request.backend = backendRequest;

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

    case {'periodic', 'periodic_ewald'}
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
             'The pair cache is an internal implementation detail of the nonperiodic matrix-free operator.']);

    otherwise
        error('thole:make_polarization_operator:UnsupportedOperatorBackend', ...
            ['Unsupported Backend="%s". Supported public backend values are: ', ...
             '"auto", "dense", "matrix_free".'], backend);
end

end

function backend = local_resolve_backend(mode, solver, backendRequest, rcut)

if ~strcmp(mode, 'nonperiodic')
    error('thole:make_polarization_operator:UnsupportedOperatorMode', ...
        'Only mode="nonperiodic" is currently supported.');
end

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
            ['Nonperiodic matrix-free operator construction currently requires a finite Rcut. ', ...
             'Use Backend="dense" for full all-pairs calculations, or provide a finite Rcut.']);
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
                ['Backend="auto" with Solver="%s" resolves to a matrix-free operator, ', ...
                 'which currently requires a finite Rcut for nonperiodic systems. ', ...
                 'Use Backend="dense" for full all-pairs calculations, or provide a finite Rcut.'], solver);
        end

        backend = 'matrix_free_paircache';

    case 'sor'
        if ~isfinite(rcut)
            error('thole:make_polarization_operator:SorRequiresFiniteCutoff', ...
                ['Backend="auto" with Solver="sor" resolves to a matrix-free row-cache operator, ', ...
                 'which currently requires a finite Rcut for nonperiodic systems. ', ...
                 'Provide a finite Rcut.']);
        end

        backend = 'matrix_free_rowcache';

    otherwise
        error('thole:make_polarization_operator:UnsupportedSolver', ...
            'Unsupported solver "%s".', solver);
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

end