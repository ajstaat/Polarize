function [mu, info] = solve_scf_gmres(problem, opOrTpol, opts)
%SOLVE_SCF_GMRES Krylov/GMRES SCF solver using operator apply().
%
% [mu, info] = thole.solve_scf_gmres(problem, op)
% [mu, info] = thole.solve_scf_gmres(problem, Tpol)
% [mu, info] = thole.solve_scf_gmres(problem, opOrTpol, opts)
%
% Solves the active-space linear equation:
%
%   (I - A*T) mu = A*Eext
%
% using MATLAB gmres and only requiring an operator action T*mu.
%
% Optional opts fields:
%   tol
%       GMRES relative tolerance, default problem.tol or 1e-10
%
%   max_iter
%       maximum GMRES iterations. For unrestarted GMRES, this is bounded by
%       the active-space dimension before calling MATLAB gmres.
%       default: problem.maxIter or 100
%
%   restart
%       restart length. Default [] means unrestarted GMRES. For larger
%       systems, a value like 30 or 50 is often safer.
%
%   verbose
%       logical, default false
%
% Output:
%   mu
%       full-system N x 3 induced dipoles
%
%   info
%       diagnostic struct

if nargin < 3 || isempty(opts)
    opts = struct();
end

validate_problem(problem);

op = local_normalize_operator(opOrTpol);

tol = local_get_field(opts, 'tol', local_get_field(problem, 'tol', 1e-10));
maxIterRequested = local_get_field(opts, 'max_iter', local_get_field(problem, 'maxIter', 100));
restart = local_get_field(opts, 'restart', []);
verbose = local_get_field(opts, 'verbose', false);

validate_options(tol, maxIterRequested, restart, verbose);

alphaVec = problem.alpha_pol_vec(:);
Evec = problem.Eext_pol_vec(:);
nVec = numel(Evec);

rhs = alphaVec .* Evec;

Afun = @(x) local_apply_linear_system(x, alphaVec, op);

if isfield(problem, 'mu0_pol_vec') && ~isempty(problem.mu0_pol_vec)
    x0 = problem.mu0_pol_vec(:);
else
    x0 = zeros(nVec, 1);
end

if numel(x0) ~= nVec
    error('thole:solve_scf_gmres:BadInitialGuess', ...
        'Initial active-space dipole vector has wrong size.');
end

% MATLAB warns if unrestarted GMRES receives maxit > size(A,1), then
% silently clamps it. Do that ourselves so tests/workflows are quiet and the
% effective value is explicit in info.
if isempty(restart)
    maxIterEffective = min(maxIterRequested, nVec);
else
    maxIterEffective = maxIterRequested;
end

if verbose
    fprintf('solve_scf_gmres:\n');
    fprintf('  nPolSites = %d\n', problem.nPolSites);
    fprintf('  nVec      = %d\n', nVec);
    fprintf('  tol       = %.3e\n', tol);

    fprintf('  max_iter  = %d', maxIterEffective);
    if maxIterEffective ~= maxIterRequested
        fprintf(' (requested %d)', maxIterRequested);
    end
    fprintf('\n');

    if isempty(restart)
        fprintf('  restart   = [] unrestarted\n');
    else
        fprintf('  restart   = %d\n', restart);
    end

    fprintf('  op.kind   = %s\n', op.kind);
    fprintf('  op.backend= %s\n', op.backend);
end

tStart = tic;

% MATLAB gmres accepts:
%   gmres(Afun, b, restart, tol, maxit, M1, M2, x0)
%
% If restart = [], gmres is unrestarted and maxit is the maximum number of
% iterations. For restarted GMRES, maxit is the maximum number of outer
% iterations.
[muPolVec, flag, gmresRelres, iter, resvec] = gmres( ...
    Afun, rhs, restart, tol, maxIterEffective, [], [], x0);

solveTime = toc(tStart);

mu = zeros(problem.nSites, 3);
mu(problem.activeSites, :) = util.unstack_xyz(muPolVec);

% Compute one explicit fixed-point residual for consistency with the other
% solvers. This uses one extra op.apply call.
resVecFinal = muPolVec - alphaVec .* (Evec + op.apply(muPolVec));
rhsScale = norm(rhs);
if rhsScale == 0
    rhsScale = 1;
end
scfRelres = norm(resVecFinal) / rhsScale;

info = struct();
info.method = 'gmres';
info.converged = (flag == 0);
info.flag = flag;

info.relres = scfRelres;
info.gmres_relres = gmresRelres;

info.iter = iter;
info.resvec = resvec;

info.tol = tol;
info.max_iter_requested = maxIterRequested;
info.max_iter = maxIterEffective;
info.restart = restart;

info.solve_time = solveTime;
info.nPolSites = problem.nPolSites;
info.nActiveVec = nVec;

info.operator_kind = op.kind;
info.operator_backend = op.backend;

if verbose
    fprintf('  flag      = %d\n', flag);
    fprintf('  gmres rel = %.3e\n', gmresRelres);
    fprintf('  scf rel   = %.3e\n', scfRelres);
    fprintf('  time      = %.6f s\n', solveTime);
end

end

% =========================================================================
% Core operator
% =========================================================================

function y = local_apply_linear_system(x, alphaVec, op)

x = x(:);
Tmu = op.apply(x);

% (I - A*T) x
y = x - alphaVec .* Tmu;

end

% =========================================================================
% Helpers
% =========================================================================

function validate_problem(problem)

if ~isstruct(problem)
    error('thole:solve_scf_gmres:BadProblem', ...
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
        error('thole:solve_scf_gmres:MissingProblemField', ...
            'problem.%s is required and missing/empty.', name);
    end
end

nVec = 3 * problem.nPolSites;

if numel(problem.Eext_pol_vec) ~= nVec || numel(problem.alpha_pol_vec) ~= nVec
    error('thole:solve_scf_gmres:BadProblemSize', ...
        'Active-space vectors must have length 3*problem.nPolSites.');
end

if numel(problem.activeSites) ~= problem.nPolSites
    error('thole:solve_scf_gmres:BadActiveSites', ...
        'numel(problem.activeSites) must equal problem.nPolSites.');
end

end

function validate_options(tol, maxIter, restart, verbose)

if ~(isnumeric(tol) && isscalar(tol) && isfinite(tol) && tol > 0)
    error('thole:solve_scf_gmres:BadTol', ...
        'tol must be a positive finite scalar.');
end

if ~(isnumeric(maxIter) && isscalar(maxIter) && maxIter >= 1 && maxIter == round(maxIter))
    error('thole:solve_scf_gmres:BadMaxIter', ...
        'max_iter must be a positive integer.');
end

if ~isempty(restart)
    if ~(isnumeric(restart) && isscalar(restart) && restart >= 1 && restart == round(restart))
        error('thole:solve_scf_gmres:BadRestart', ...
            'restart must be empty or a positive integer.');
    end
end

if ~(islogical(verbose) && isscalar(verbose))
    error('thole:solve_scf_gmres:BadVerbose', ...
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
    error('thole:solve_scf_gmres:BadOperator', ...
        'Second input must be a Tpol matrix or operator struct.');
end

if ~isfield(opOrTpol, 'apply') || ~isa(opOrTpol.apply, 'function_handle')
    backend = '<unknown>';
    if isfield(opOrTpol, 'backend')
        backend = char(string(opOrTpol.backend));
    end

    error('thole:solve_scf_gmres:OperatorMissingApply', ...
        ['GMRES SCF requires an operator with op.apply(muVec). ', ...
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

function value = local_get_field(s, name, defaultValue)

if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end

end