function [mu, info] = solve_scf_direct(problem, opOrTpol, opts)
%SOLVE_SCF_DIRECT Solve linear Thole SCF problem by dense direct solve.
%
% [mu, info] = thole.solve_scf_direct(problem, Tpol)
% [mu, info] = thole.solve_scf_direct(problem, op)
% [mu, info] = thole.solve_scf_direct(problem, opOrTpol, opts)
%
% Solves the active-space equation:
%
%   (I - A*T) mu = A*Eext
%
% where:
%
%   A       = diag(problem.alpha_pol_vec)
%   T       = active-space dipole interaction matrix
%   Eext    = problem.Eext_pol_vec
%   mu      = active-space induced dipole vector
%
% The returned mu is expanded to the full system size as N x 3, with zeros
% on non-polarizable/non-active sites.
%
% Optional opts fields:
%   compute_residual
%       logical, default true.
%       If true, compute thole.compute_active_space_relres after solving.
%       This is useful for validation but can be expensive for large dense
%       matrices because it requires an additional Tpol * mu multiply.
%
% Output info fields:
%   method
%   nPolSites
%   nActiveVec
%   setup_time
%   solve_time
%   residual_time
%   total_time
%   compute_residual
%   relres
%   converged
%   tol
%
% Notes
%   This routine intentionally does not compute cond(M) or rcond(M).
%   Those diagnostics are expensive for large dense active-space matrices.
%
%   It also intentionally avoids forming A = diag(alphaVec). Instead:
%
%       A*T    = alphaVec .* T
%       A*Eext = alphaVec .* Eext
%
%   where alphaVec is a column vector and implicit expansion row-scales T.

if nargin < 3 || isempty(opts)
    opts = struct();
end

tTotal = tic;

validate_problem(problem);

computeResidual = local_get_field(opts, 'compute_residual', true);

if ~(islogical(computeResidual) && isscalar(computeResidual))
    error('thole:solve_scf_direct:BadComputeResidual', ...
        'opts.compute_residual must be a logical scalar.');
end

Tpol = local_extract_dense_Tpol(opOrTpol);

nVec = numel(problem.Eext_pol_vec);

if ~isequal(size(Tpol), [nVec nVec])
    error('thole:solve_scf_direct:BadTpolSize', ...
        'Tpol must be %d x %d for this active-space problem.', nVec, nVec);
end

alphaVec = problem.alpha_pol_vec(:);
Evec = problem.Eext_pol_vec(:);

% -------------------------------------------------------------------------
% Matrix/rhs setup
% -------------------------------------------------------------------------

tSetup = tic;

% Row-scaled interaction matrix:
%
%   diag(alphaVec) * Tpol
%
% without explicitly forming diag(alphaVec).
AT = alphaVec .* Tpol;

M = eye(nVec) - AT;
rhs = alphaVec .* Evec;

setupTime = toc(tSetup);

% -------------------------------------------------------------------------
% Direct dense linear solve
% -------------------------------------------------------------------------

tSolve = tic;
muPolVec = M \ rhs;
solveTime = toc(tSolve);

mu = zeros(problem.nSites, 3);
mu(problem.activeSites, :) = util.unstack_xyz(muPolVec);

% -------------------------------------------------------------------------
% Optional residual diagnostic
% -------------------------------------------------------------------------

relres = NaN;
residualTime = 0.0;

if computeResidual
    tRes = tic;
    relres = thole.compute_active_space_relres(problem, Tpol, mu);
    residualTime = toc(tRes);
end

totalTime = toc(tTotal);

info = struct();
info.method = 'direct';
info.nPolSites = problem.nPolSites;
info.nActiveVec = nVec;

info.setup_time = setupTime;
info.solve_time = solveTime;
info.residual_time = residualTime;
info.total_time = totalTime;

info.compute_residual = computeResidual;
info.relres = relres;

if computeResidual
    info.converged = relres <= problem.tol;
else
    info.converged = [];
end

info.tol = problem.tol;

end

% =========================================================================
% Helpers
% =========================================================================

function validate_problem(problem)

if ~isstruct(problem)
    error('thole:solve_scf_direct:BadProblem', ...
        'problem must be a struct from thole.prepare_scf_problem.');
end

required = {
    'activeSites'
    'nPolSites'
    'nSites'
    'Eext_pol_vec'
    'alpha_pol_vec'
    'tol'
};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(problem, name) || isempty(problem.(name))
        error('thole:solve_scf_direct:MissingProblemField', ...
            'problem.%s is required and missing/empty.', name);
    end
end

if ~(isnumeric(problem.nPolSites) && isscalar(problem.nPolSites) && problem.nPolSites >= 0)
    error('thole:solve_scf_direct:BadNPolSites', ...
        'problem.nPolSites must be a nonnegative scalar.');
end

if ~(isnumeric(problem.nSites) && isscalar(problem.nSites) && problem.nSites >= problem.nPolSites)
    error('thole:solve_scf_direct:BadNSites', ...
        'problem.nSites must be a scalar >= problem.nPolSites.');
end

if numel(problem.activeSites) ~= problem.nPolSites
    error('thole:solve_scf_direct:BadActiveSites', ...
        'numel(problem.activeSites) must equal problem.nPolSites.');
end

nVec = 3 * problem.nPolSites;

if numel(problem.Eext_pol_vec) ~= nVec
    error('thole:solve_scf_direct:BadEextSize', ...
        'problem.Eext_pol_vec must have length 3*problem.nPolSites.');
end

if numel(problem.alpha_pol_vec) ~= nVec
    error('thole:solve_scf_direct:BadAlphaSize', ...
        'problem.alpha_pol_vec must have length 3*problem.nPolSites.');
end

if any(problem.alpha_pol_vec(:) < 0)
    error('thole:solve_scf_direct:NegativeAlpha', ...
        'problem.alpha_pol_vec must be nonnegative.');
end

end

function Tpol = local_extract_dense_Tpol(opOrTpol)

if isnumeric(opOrTpol)
    Tpol = opOrTpol;
    return;
end

if ~isstruct(opOrTpol)
    error('thole:solve_scf_direct:BadOperator', ...
        'Second input must be a dense Tpol matrix or an operator struct.');
end

if ~isfield(opOrTpol, 'kind') || ~strcmp(opOrTpol.kind, 'dense_matrix') || ...
        ~isfield(opOrTpol, 'Tpol') || isempty(opOrTpol.Tpol)

    backend = '<unknown>';
    if isfield(opOrTpol, 'backend') && ~isempty(opOrTpol.backend)
        backend = char(string(opOrTpol.backend));
    end

    error('thole:solve_scf_direct:RequiresDenseOperator', ...
        ['Direct SCF requires a dense polarization operator with op.Tpol. ', ...
         'Requested backend "%s" is not dense; use Backend="dense" or choose ', ...
         'an iterative solver that supports matrix-free operators.'], backend);
end

Tpol = opOrTpol.Tpol;

if ~isnumeric(Tpol) || ndims(Tpol) ~= 2 || size(Tpol,1) ~= size(Tpol,2)
    error('thole:solve_scf_direct:BadDenseOperatorMatrix', ...
        'op.Tpol must be a numeric square matrix.');
end

end

function value = local_get_field(s, name, defaultValue)

if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end

end