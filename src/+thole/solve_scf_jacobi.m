function [mu, info] = solve_scf_jacobi(problem, opOrTpol, opts)
%SOLVE_SCF_JACOBI Jacobi/fixed-point SCF solver using operator apply().
%
% [mu, info] = thole.solve_scf_jacobi(problem, op)
% [mu, info] = thole.solve_scf_jacobi(problem, Tpol)
% [mu, info] = thole.solve_scf_jacobi(problem, opOrTpol, opts)
%
% Fixed-point equation:
%
%   mu = A * (Eext + T*mu)
%
% with linear mixing:
%
%   mu_{n+1} = (1-mix) * mu_n + mix * A*(Eext + T*mu_n)
%
% This solver requires only op.apply(muVec), so it can use dense or
% matrix-free operators.
%
% Options:
%   tol
%       convergence tolerance, default problem.tol or 1e-10
%
%   max_iter
%       maximum iterations, default problem.maxIter or 500
%
%   mixing
%       linear mixing parameter in (0,1], default problem.mixing or 0.5
%
%   stop_metric
%       'relres'  : stop on relative SCF residual
%       'max_dmu' : stop on max per-site dipole update magnitude
%       default: problem.stop_metric/problem.stopMetric or 'relres'
%
%   verbose
%       logical, default false

if nargin < 3 || isempty(opts)
    opts = struct();
end

validate_problem(problem);

op = local_normalize_operator(opOrTpol);

tol = local_get_field(opts, 'tol', local_get_field(problem, 'tol', 1e-10));
maxIter = local_get_field(opts, 'max_iter', local_get_field(problem, 'maxIter', 500));
mixing = local_get_field(opts, 'mixing', local_get_field(problem, 'mixing', 0.5));
stopMetric = local_normalize_stop_metric(local_get_field(opts, 'stop_metric', 'relres'));
verbose = local_get_field(opts, 'verbose', false);

validate_options(tol, maxIter, mixing, verbose);

alphaVec = problem.alpha_pol_vec(:);
Evec = problem.Eext_pol_vec(:);
nVec = numel(Evec);

if isfield(problem, 'mu0_pol_vec') && ~isempty(problem.mu0_pol_vec)
    muVec = problem.mu0_pol_vec(:);
else
    muVec = zeros(nVec, 1);
end

if numel(muVec) ~= nVec
    error('thole:solve_scf_jacobi:BadInitialGuess', ...
        'Initial active-space dipole vector has wrong size.');
end

relresHistory = NaN(maxIter, 1);
deltaHistory = NaN(maxIter, 1);
maxDmuHistory = NaN(maxIter, 1);
stopHistory = NaN(maxIter, 1);

rhsScale = norm(alphaVec .* Evec);
if rhsScale == 0
    rhsScale = 1;
end

tStart = tic;

converged = false;
iter = 0;
relres = Inf;
delta = Inf;
maxDmu = Inf;
stopValue = Inf;

for k = 1:maxIter
    iter = k;

    Tmu = op.apply(muVec);
    muFixed = alphaVec .* (Evec + Tmu);

    muNew = (1 - mixing) * muVec + mixing * muFixed;

    % Residual of updated iterate.
    TmuNew = op.apply(muNew);
    resVec = muNew - alphaVec .* (Evec + TmuNew);

    relres = norm(resVec) / rhsScale;
    delta = norm(muNew - muVec) / max(norm(muNew), eps);

    dmuVec = muNew - muVec;
    dmuMat = util.unstack_xyz(dmuVec);
    maxDmu = max(vecnorm(dmuMat, 2, 2));

    switch stopMetric
        case 'relres'
            stopValue = relres;
        case 'max_dmu'
            stopValue = maxDmu;
        otherwise
            error('thole:solve_scf_jacobi:BadStopMetric', ...
                'Unsupported stop metric "%s".', stopMetric);
    end

    relresHistory(k) = relres;
    deltaHistory(k) = delta;
    maxDmuHistory(k) = maxDmu;
    stopHistory(k) = stopValue;

    muVec = muNew;

    if verbose
        fprintf('  Jacobi iter %4d: relres = %.3e, max_dmu = %.3e, delta = %.3e\n', ...
            k, relres, maxDmu, delta);
    end

    if stopValue <= tol
        converged = true;
        break;
    end
end

solveTime = toc(tStart);

mu = zeros(problem.nSites, 3);
mu(problem.activeSites, :) = util.unstack_xyz(muVec);

info = struct();
info.method = 'jacobi';
info.converged = converged;
info.iterations = iter;

info.relres = relres;
info.delta = delta;
info.max_dmu = maxDmu;

info.stop_metric = stopMetric;
info.stop_value = stopValue;

info.tol = tol;
info.max_iter = maxIter;
info.mixing = mixing;
info.solve_time = solveTime;
info.nPolSites = problem.nPolSites;
info.nActiveVec = nVec;

info.relres_history = relresHistory(1:iter);
info.delta_history = deltaHistory(1:iter);
info.max_dmu_history = maxDmuHistory(1:iter);
info.stop_history = stopHistory(1:iter);

info.operator_kind = op.kind;
info.operator_backend = op.backend;

end

% =========================================================================
% Helpers
% =========================================================================

function validate_problem(problem)

if ~isstruct(problem)
    error('thole:solve_scf_jacobi:BadProblem', ...
        'problem must be a struct from thole.prepare_scf_problem.');
end

required = {
    'activeSites'
    'nPolSites'
    'nSites'
    'Eext_pol_vec'
    'alpha_pol_vec'
};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(problem, name) || isempty(problem.(name))
        error('thole:solve_scf_jacobi:MissingProblemField', ...
            'problem.%s is required and missing/empty.', name);
    end
end

nVec = 3 * problem.nPolSites;

if numel(problem.Eext_pol_vec) ~= nVec || numel(problem.alpha_pol_vec) ~= nVec
    error('thole:solve_scf_jacobi:BadProblemSize', ...
        'Active-space vectors must have length 3*problem.nPolSites.');
end

if numel(problem.activeSites) ~= problem.nPolSites
    error('thole:solve_scf_jacobi:BadActiveSites', ...
        'numel(problem.activeSites) must equal problem.nPolSites.');
end

end

function validate_options(tol, maxIter, mixing, verbose)

if ~(isnumeric(tol) && isscalar(tol) && isfinite(tol) && tol > 0)
    error('thole:solve_scf_jacobi:BadTol', ...
        'tol must be a positive finite scalar.');
end

if ~(isnumeric(maxIter) && isscalar(maxIter) && maxIter >= 1 && maxIter == round(maxIter))
    error('thole:solve_scf_jacobi:BadMaxIter', ...
        'max_iter must be a positive integer.');
end

if ~(isnumeric(mixing) && isscalar(mixing) && isfinite(mixing) && mixing > 0 && mixing <= 1)
    error('thole:solve_scf_jacobi:BadMixing', ...
        'mixing must be in (0, 1].');
end

if ~(islogical(verbose) && isscalar(verbose))
    error('thole:solve_scf_jacobi:BadVerbose', ...
        'verbose must be a logical scalar.');
end

end

function op = local_normalize_operator(opOrTpol)

if isnumeric(opOrTpol)
    Tpol = opOrTpol;

    op = struct();
    op.mode = 'unknown';
    op.kind = 'dense_matrix';
    op.backend = 'raw_dense_matrix';
    op.Tpol = Tpol;
    op.apply = @(muVec) Tpol * muVec;
    return;
end

if ~isstruct(opOrTpol)
    error('thole:solve_scf_jacobi:BadOperator', ...
        'Second input must be a Tpol matrix or operator struct.');
end

if ~isfield(opOrTpol, 'apply') || ~isa(opOrTpol.apply, 'function_handle')
    backend = '<unknown>';
    if isfield(opOrTpol, 'backend')
        backend = char(string(opOrTpol.backend));
    end

    error('thole:solve_scf_jacobi:OperatorMissingApply', ...
        ['Jacobi SCF requires an operator with op.apply(muVec). ', ...
         'Requested backend "%s" does not provide apply().'], backend);
end

op = opOrTpol;

if ~isfield(op, 'backend')
    op.backend = '<unknown>';
end

if ~isfield(op, 'kind')
    op.kind = '<unknown>';
end

end

function stopMetric = local_normalize_stop_metric(x)

stopMetric = lower(strtrim(char(string(x))));

switch stopMetric
    case {'relres', 'relative_residual', 'residual'}
        stopMetric = 'relres';

    case {'max_dmu', 'max_dipole_update', 'dmu', 'max_update'}
        stopMetric = 'max_dmu';

    otherwise
        error('thole:solve_scf_jacobi:BadStopMetric', ...
            'stop_metric must be ''relres'' or ''max_dmu''.');
end

end

function value = local_get_field(s, name, defaultValue)

if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end

end