function op = make_polarization_operator(sys, problem, varargin)
%MAKE_POLARIZATION_OPERATOR Build a polarization operator object.
%
% op = thole.make_polarization_operator(sys, problem, Name, Value)
%
% The operator object standardizes how solvers access the dipole-dipole
% interaction operator T. Dense, matrix-free, periodic, and future P3M
% backends should all present a consistent interface.
%
% Currently supported:
%
%   Mode    = 'nonperiodic'
%   Backend = 'dense'
%
% Options
%   'Mode'       'nonperiodic' | 'periodic_ewald' | 'periodic_p3m'
%                default: inferred from sys.is_periodic, otherwise nonperiodic
%
%   'Backend'    'dense' | future backends
%                default: 'dense'
%
%   'UseThole'   logical, default true
%   'Softening'  scalar, default 0
%   'Rcut'       scalar cutoff in bohr, default Inf
%   'Verbose'    logical, default false
%
% Output op fields
%   .mode
%   .backend
%   .kind
%   .nPolSites
%   .size
%   .apply
%   .info
%
% Dense backend additionally contains:
%   .Tpol

p = inputParser;
addRequired(p, 'sys', @isstruct);
addRequired(p, 'problem', @isstruct);
addParameter(p, 'Mode', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Backend', 'dense', @(x) ischar(x) || isstring(x));
addParameter(p, 'UseThole', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Softening', 0.0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'Rcut', Inf, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
parse(p, sys, problem, varargin{:});

opt = p.Results;

mode = lower(strtrim(char(string(opt.Mode))));
backend = lower(strtrim(char(string(opt.Backend))));

if isempty(mode)
    mode = local_infer_mode(sys);
end

validate_problem(problem);

switch mode
    case {'nonperiodic', 'finite', 'cluster'}
        mode = 'nonperiodic';

    case {'periodic', 'periodic_ewald'}
        mode = 'periodic_ewald';

    case {'periodic_p3m', 'p3m'}
        mode = 'periodic_p3m';

    otherwise
        error('thole:make_polarization_operator:UnsupportedMode', ...
            ['Unsupported polarization-operator mode "%s". Supported modes are: ', ...
             '"nonperiodic", "periodic_ewald", and future "periodic_p3m".'], mode);
end

switch backend
    case {'dense', 'dense_matrix'}
        backend = 'dense';

    case {'pair_cache', 'matrix_free', 'ewald_direct', 'p3m'}
        % Recognized names, but not implemented in Operator-1.
        error('thole:make_polarization_operator:UnsupportedOperatorBackend', ...
            ['Unsupported polarization operator mode/backend combination: ', ...
             'mode="%s", backend="%s". Currently supported: ', ...
             'mode="nonperiodic", backend="dense".'], mode, backend);

    otherwise
        error('thole:make_polarization_operator:UnsupportedOperatorBackend', ...
            ['Unsupported polarization-operator backend "%s". Currently supported: ', ...
             'backend="dense" for mode="nonperiodic".'], backend);
end

if strcmp(mode, 'nonperiodic') && strcmp(backend, 'dense')
    scfParams = struct();
    scfParams.use_thole = opt.UseThole;
    scfParams.softening = opt.Softening;
    scfParams.rcut = opt.Rcut;
    scfParams.verbose = opt.Verbose;

    [Tpol, info] = thole.assemble_nonperiodic_interaction_matrix(sys, problem, scfParams);

    op = struct();
    op.mode = 'nonperiodic';
    op.backend = 'dense';
    op.kind = 'dense_matrix';
    op.nPolSites = problem.nPolSites;
    op.size = size(Tpol);
    op.Tpol = Tpol;
    op.apply = @(muVec) Tpol * muVec;
    op.info = info;

    op.params = struct();
    op.params.use_thole = opt.UseThole;
    op.params.softening = opt.Softening;
    op.params.rcut = opt.Rcut;

    return;
end

error('thole:make_polarization_operator:UnsupportedOperatorBackend', ...
    ['Unsupported polarization operator mode/backend combination: ', ...
     'mode="%s", backend="%s". Currently supported: ', ...
     'mode="nonperiodic", backend="dense".'], mode, backend);

end

% =========================================================================
% Helpers
% =========================================================================

function mode = local_infer_mode(sys)

if isfield(sys, 'is_periodic') && ~isempty(sys.is_periodic) && logical(sys.is_periodic)
    mode = 'periodic_ewald';
else
    mode = 'nonperiodic';
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