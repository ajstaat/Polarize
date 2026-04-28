function [mu, info] = solve_scf_sor(problem, op, opts)
%SOLVE_SCF_SOR Site-block SOR/Gauss-Seidel SCF solver.
%
% [mu, info] = thole.solve_scf_sor(problem, op)
% [mu, info] = thole.solve_scf_sor(problem, op, opts)
%
% Solves:
%
%   mu = A * (Eext + T*mu)
%
% using active-site block Gauss-Seidel / SOR updates:
%
%   mu_i <- (1-omega)*mu_i + omega*alpha_i*(E_i + sum_j T_ij mu_j)
%
% Requirements:
%   op.kind = 'matrix_free'
%   op.capabilities.row_update = true
%   op.apply_row(rowLocal, muVec)
%   op.apply(muVec)
%
% Options:
%   tol       default problem.tol or 1e-10
%   max_iter  default problem.maxIter or 500
%   omega     default problem.omega or 1.0
%   verbose   default false

if nargin < 3 || isempty(opts)
    opts = struct();
end

validate_problem(problem);
validate_operator(op);

tol = local_get_field(opts, 'tol', local_get_field(problem, 'tol', 1e-10));
maxIter = local_get_field(opts, 'max_iter', local_get_field(problem, 'maxIter', 500));
omega = local_get_field(opts, 'omega', local_get_field(problem, 'omega', 1.0));
verbose = local_get_field(opts, 'verbose', false);

validate_options(tol, maxIter, omega, verbose);

nPol = problem.nPolSites;
nVec = 3*nPol;

alphaVec = problem.alpha_pol_vec(:);
Evec = problem.Eext_pol_vec(:);

if isfield(problem, 'mu0_pol_vec') && ~isempty(problem.mu0_pol_vec)
    muVec = problem.mu0_pol_vec(:);
else
    muVec = zeros(nVec, 1);
end

if numel(muVec) ~= nVec
    error('thole:solve_scf_sor:BadInitialGuess', ...
        'Initial active-space dipole vector has wrong size.');
end

rhsScale = norm(alphaVec .* Evec);
if rhsScale == 0
    rhsScale = 1;
end

relresHistory = NaN(maxIter, 1);
deltaHistory = NaN(maxIter, 1);

tStart = tic;

converged = false;
relres = Inf;
iter = 0;

for k = 1:maxIter
    iter = k;

    muOld = muVec;

    for a = 1:nPol
        block = local_block_indices(a);

        TiMu = op.apply_row(a, muVec);

        muFixed = alphaVec(block) .* (Evec(block) + TiMu);

        muVec(block) = (1 - omega) .* muVec(block) + omega .* muFixed;
    end

    Tmu = op.apply(muVec);
    resVec = muVec - alphaVec .* (Evec + Tmu);

    relres = norm(resVec) / rhsScale;
    delta = norm(muVec - muOld) / max(norm(muVec), eps);

    relresHistory(k) = relres;
    deltaHistory(k) = delta;

    if verbose
        fprintf('  SOR iter %4d: relres = %.3e, delta = %.3e\n', k, relres, delta);
    end

    if relres <= tol
        converged = true;
        break;
    end
end

solveTime = toc(tStart);

mu = zeros(problem.nSites, 3);
mu(problem.activeSites, :) = util.unstack_xyz(muVec);

info = struct();
info.method = 'sor';
info.converged = converged;
info.iterations = iter;
info.relres = relres;
info.tol = tol;
info.max_iter = maxIter;
info.omega = omega;
info.solve_time = solveTime;
info.nPolSites = nPol;
info.nActiveVec = nVec;
info.relres_history = relresHistory(1:iter);
info.delta_history = deltaHistory(1:iter);
info.operator_kind = op.kind;
info.operator_backend = op.backend;

end

% =========================================================================
% Helpers
% =========================================================================

function validate_problem(problem)

if ~isstruct(problem)
    error('thole:solve_scf_sor:BadProblem', ...
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
        error('thole:solve_scf_sor:MissingProblemField', ...
            'problem.%s is required and missing/empty.', name);
    end
end

nVec = 3 * problem.nPolSites;

if numel(problem.Eext_pol_vec) ~= nVec || numel(problem.alpha_pol_vec) ~= nVec
    error('thole:solve_scf_sor:BadProblemSize', ...
        'Active-space vectors must have length 3*problem.nPolSites.');
end

if numel(problem.activeSites) ~= problem.nPolSites
    error('thole:solve_scf_sor:BadActiveSites', ...
        'numel(problem.activeSites) must equal problem.nPolSites.');
end

end

function validate_operator(op)

if ~isstruct(op)
    error('thole:solve_scf_sor:BadOperator', ...
        'SOR requires an operator struct.');
end

if ~isfield(op, 'kind') || ~strcmp(op.kind, 'matrix_free')
    kind = '<unknown>';
    if isfield(op, 'kind')
        kind = char(string(op.kind));
    end

    error('thole:solve_scf_sor:RequiresMatrixFreeOperator', ...
        ['SOR requires a matrix-free row-update operator. ', ...
         'Requested op.kind="%s".'], kind);
end

if ~isfield(op, 'capabilities') || ~isfield(op.capabilities, 'row_update') || ...
        ~op.capabilities.row_update
    error('thole:solve_scf_sor:RequiresRowUpdateCapability', ...
        ['SOR requires op.capabilities.row_update = true. ', ...
         'Construct the operator with Solver="sor", Backend="auto" or "matrix_free".']);
end

if ~isfield(op, 'apply_row') || ~isa(op.apply_row, 'function_handle')
    error('thole:solve_scf_sor:MissingApplyRow', ...
        'SOR requires op.apply_row(rowLocal, muVec).');
end

if ~isfield(op, 'apply') || ~isa(op.apply, 'function_handle')
    error('thole:solve_scf_sor:MissingApply', ...
        'SOR requires op.apply(muVec) for residual diagnostics.');
end

end

function validate_options(tol, maxIter, omega, verbose)

if ~(isnumeric(tol) && isscalar(tol) && isfinite(tol) && tol > 0)
    error('thole:solve_scf_sor:BadTol', ...
        'tol must be a positive finite scalar.');
end

if ~(isnumeric(maxIter) && isscalar(maxIter) && maxIter >= 1 && maxIter == round(maxIter))
    error('thole:solve_scf_sor:BadMaxIter', ...
        'max_iter must be a positive integer.');
end

if ~(isnumeric(omega) && isscalar(omega) && isfinite(omega) && omega > 0 && omega < 2)
    error('thole:solve_scf_sor:BadOmega', ...
        'omega must be in the interval (0, 2).');
end

if ~(islogical(verbose) && isscalar(verbose))
    error('thole:solve_scf_sor:BadVerbose', ...
        'verbose must be a logical scalar.');
end

end

function idx = local_block_indices(k)

idx = (3*(k-1)+1):(3*k);

end

function value = local_get_field(s, name, defaultValue)

if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end

end