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
% Fast path:
%   If op.row_cache contains raw row-cache arrays:
%
%       row_ptr, col_idx, dr, thole_f3, thole_f5,
%       inv_r3_bare, inv_r5_bare
%
%   this solver performs the old-style matrix-free SOR contraction directly
%   on those arrays. This avoids one function-handle call per active site
%   per iteration and avoids the slower tensor-coefficient SOR path.
%
% Options:
%   tol
%       convergence tolerance, default problem.tol or 1e-10
%
%   max_iter
%       maximum iterations, default problem.maxIter or 500
%
%   omega
%       SOR relaxation parameter in (0,2), default problem.omega or 1.0
%
%   stop_metric
%       'relres'  : stop on relative SCF residual
%       'max_dmu' : stop on max per-site dipole update magnitude
%       default: 'relres'
%
%   residual_every
%       cadence for expensive full residual evaluations when
%       stop_metric='max_dmu'.
%
%       1      : compute residual every iteration
%       25     : compute residual on iter 1, 25, 50, ...
%       0/Inf  : only compute final residual after convergence/end
%
%       If stop_metric='relres', residuals are computed every iteration
%       regardless of residual_every.
%
%   verbose
%       logical, default false

if nargin < 3 || isempty(opts)
    opts = struct();
end

validate_problem(problem);
validate_operator(op);

tol = local_get_field(opts, 'tol', local_get_field(problem, 'tol', 1e-10));
maxIter = local_get_field(opts, 'max_iter', local_get_field(problem, 'maxIter', 500));
omega = local_get_field(opts, 'omega', local_get_field(problem, 'omega', 1.0));

% Keep solver default clean and explicit. Legacy problem.stopMetric is not
% inherited implicitly; callers should pass opts.stop_metric when they want
% old-style max_dmu stopping.
stopMetric = local_normalize_stop_metric(local_get_field(opts, 'stop_metric', 'relres'));

residualEvery = local_get_field(opts, 'residual_every', ...
    local_get_field(problem, 'residual_every', ...
    local_get_field(problem, 'residualEvery', 1)));

verbose = local_get_field(opts, 'verbose', false);

validate_options(tol, maxIter, omega, residualEvery, verbose);

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

useRawRowCacheFastPath = local_has_raw_rowcache_fast_path(op);

rhsScale = norm(alphaVec .* Evec);
if rhsScale == 0
    rhsScale = 1;
end

relresHistory = NaN(maxIter, 1);
deltaHistory = NaN(maxIter, 1);
maxDmuHistory = NaN(maxIter, 1);
stopHistory = NaN(maxIter, 1);
residualComputedHistory = false(maxIter, 1);

tStart = tic;

converged = false;
relres = Inf;
delta = Inf;
maxDmu = Inf;
stopValue = Inf;
iter = 0;

for k = 1:maxIter
    iter = k;

    muOldVec = muVec;

    % ---------------------------------------------------------------------
    % Gauss-Seidel / SOR sweep.
    %
    % Preferred production path is the raw row-cache contraction, matching
    % the old matrix-free SOR implementation.
    % ---------------------------------------------------------------------

    if useRawRowCacheFastPath
        muVec = local_sor_sweep_rowcache_raw( ...
            muVec, alphaVec, Evec, omega, op.row_cache);
    else
        muVec = local_sor_sweep_generic( ...
            muVec, alphaVec, Evec, omega, op, nPol);
    end

    % ---------------------------------------------------------------------
    % Cheap update diagnostics.
    % ---------------------------------------------------------------------

    dmuVec = muVec - muOldVec;
    dmuMat = util.unstack_xyz(dmuVec);

    maxDmu = max(vecnorm(dmuMat, 2, 2));
    delta = norm(dmuVec) / max(norm(muVec), eps);

    % ---------------------------------------------------------------------
    % Expensive residual diagnostic.
    %
    % If relres is the stop metric, this must be computed every iteration.
    % If max_dmu is the stop metric, compute it only at the requested
    % cadence for diagnostics, plus once after the loop for final reporting.
    % ---------------------------------------------------------------------

    needResidual = false;

    if strcmp(stopMetric, 'relres')
        needResidual = true;
    elseif local_should_compute_periodic_residual(k, residualEvery)
        needResidual = true;
    end

    if needResidual
        if useRawRowCacheFastPath
            relres = local_relres_rowcache_raw( ...
                muVec, alphaVec, Evec, op.row_cache, rhsScale);
        else
            Tmu = op.apply(muVec);
            resVec = muVec - alphaVec .* (Evec + Tmu);
            relres = norm(resVec) / rhsScale;
        end

        residualComputedHistory(k) = true;
    end

    switch stopMetric
        case 'relres'
            stopValue = relres;

        case 'max_dmu'
            stopValue = maxDmu;

        otherwise
            error('thole:solve_scf_sor:BadStopMetric', ...
                'Unsupported stop metric "%s".', stopMetric);
    end

    relresHistory(k) = relres;
    deltaHistory(k) = delta;
    maxDmuHistory(k) = maxDmu;
    stopHistory(k) = stopValue;

    if verbose
        if residualComputedHistory(k)
            fprintf('  SOR iter %4d: relres = %.3e, max_dmu = %.3e, delta = %.3e\n', ...
                k, relres, maxDmu, delta);
        else
            fprintf('  SOR iter %4d: relres = skipped, max_dmu = %.3e, delta = %.3e\n', ...
                k, maxDmu, delta);
        end
    end

    if stopValue <= tol
        converged = true;
        break;
    end
end

% Always compute a final residual for honest reporting and energy sanity.
tFinalResidual = tic;

if useRawRowCacheFastPath
    relresFinal = local_relres_rowcache_raw( ...
        muVec, alphaVec, Evec, op.row_cache, rhsScale);
else
    TmuFinal = op.apply(muVec);
    resVecFinal = muVec - alphaVec .* (Evec + TmuFinal);
    relresFinal = norm(resVecFinal) / rhsScale;
end

finalResidualTime = toc(tFinalResidual);

relres = relresFinal;

if iter > 0
    relresHistory(iter) = relresFinal;
    residualComputedHistory(iter) = true;

    if strcmp(stopMetric, 'relres')
        stopValue = relresFinal;
        stopHistory(iter) = stopValue;
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
info.delta = delta;
info.max_dmu = maxDmu;

info.stop_metric = stopMetric;
info.stop_value = stopValue;

info.tol = tol;
info.max_iter = maxIter;
info.omega = omega;
info.residual_every = residualEvery;
info.final_residual_time = finalResidualTime;

info.solve_time = solveTime;
info.nPolSites = nPol;
info.nActiveVec = nVec;

info.relres_history = relresHistory(1:iter);
info.delta_history = deltaHistory(1:iter);
info.max_dmu_history = maxDmuHistory(1:iter);
info.stop_history = stopHistory(1:iter);
info.residual_computed_history = residualComputedHistory(1:iter);

info.operator_kind = op.kind;
info.operator_backend = op.backend;
info.used_rowcache_fast_path = useRawRowCacheFastPath;
info.rowcache_fast_path_type = local_fast_path_name(useRawRowCacheFastPath);

end

% =========================================================================
% SOR sweep implementations
% =========================================================================

function muVec = local_sor_sweep_generic(muVec, alphaVec, Evec, omega, op, nPol)

for a = 1:nPol
    block = local_block_indices(a);

    TiMu = op.apply_row(a, muVec);

    muFixed = alphaVec(block) .* (Evec(block) + TiMu);

    muVec(block) = (1 - omega) .* muVec(block) + omega .* muFixed;
end

end

function muVec = local_sor_sweep_rowcache_raw(muVec, alphaVec, Evec, omega, rowCache)

rowPtr = rowCache.row_ptr(:);
colIdx = rowCache.col_idx(:);
drAll = rowCache.dr;

f3All = rowCache.thole_f3(:);
f5All = rowCache.thole_f5(:);

invR3All = rowCache.inv_r3_bare(:);
invR5All = rowCache.inv_r5_bare(:);

nPol = numel(rowPtr) - 1;

muPol = util.unstack_xyz(muVec);
muOldPol = muPol;

EextPol = util.unstack_xyz(Evec);
alphaPol = alphaVec(1:3:end);

oneMinusOmega = 1 - omega;

for i = 1:nPol
    k0 = rowPtr(i);
    k1 = rowPtr(i + 1) - 1;

    ELoc = EextPol(i, :);

    if k1 >= k0
        idx = k0:k1;
        cols = colIdx(idx);

        muNbr = muPol(cols, :);
        dr = drAll(idx, :);

        f3 = f3All(idx);
        f5 = f5All(idx);

        invR3 = invR3All(idx);
        invR5 = invR5All(idx);

        muDotR = sum(muNbr .* dr, 2);

        coeff1 = 3 .* (f5 .* muDotR .* invR5);
        coeff2 = f3 .* invR3;

        contrib = coeff1 .* dr - coeff2 .* muNbr;

        ELoc = ELoc + sum(contrib, 1);
    end

    muGS = alphaPol(i) .* ELoc;

    muPol(i, :) = oneMinusOmega .* muOldPol(i, :) + omega .* muGS;
end

muVec = util.stack_xyz(muPol);

end

% =========================================================================
% Residual implementations
% =========================================================================

function relres = local_relres_rowcache_raw(muVec, alphaVec, Evec, rowCache, rhsScale)

rowPtr = rowCache.row_ptr(:);
colIdx = rowCache.col_idx(:);
drAll = rowCache.dr;

f3All = rowCache.thole_f3(:);
f5All = rowCache.thole_f5(:);

invR3All = rowCache.inv_r3_bare(:);
invR5All = rowCache.inv_r5_bare(:);

nPol = numel(rowPtr) - 1;

muPol = util.unstack_xyz(muVec);
EextPol = util.unstack_xyz(Evec);
alphaPol = alphaVec(1:3:end);

EdipPol = zeros(nPol, 3);

for i = 1:nPol
    k0 = rowPtr(i);
    k1 = rowPtr(i + 1) - 1;

    if k1 < k0
        continue;
    end

    idx = k0:k1;
    cols = colIdx(idx);

    muNbr = muPol(cols, :);
    dr = drAll(idx, :);

    f3 = f3All(idx);
    f5 = f5All(idx);

    invR3 = invR3All(idx);
    invR5 = invR5All(idx);

    muDotR = sum(muNbr .* dr, 2);

    coeff1 = 3 .* (f5 .* muDotR .* invR5);
    coeff2 = f3 .* invR3;

    contrib = coeff1 .* dr - coeff2 .* muNbr;

    EdipPol(i, :) = sum(contrib, 1);
end

rhsPol = alphaPol .* (EextPol + EdipPol);
rPol = muPol - rhsPol;

relres = norm(rPol, 'fro') / rhsScale;

end

% =========================================================================
% Validation / option helpers
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

function validate_options(tol, maxIter, omega, residualEvery, verbose)

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

if ~(isnumeric(residualEvery) && isscalar(residualEvery) && residualEvery >= 0)
    error('thole:solve_scf_sor:BadResidualEvery', ...
        'residual_every must be a nonnegative scalar or Inf.');
end

if isfinite(residualEvery) && residualEvery ~= round(residualEvery)
    error('thole:solve_scf_sor:BadResidualEvery', ...
        'finite residual_every must be an integer.');
end

if ~(islogical(verbose) && isscalar(verbose))
    error('thole:solve_scf_sor:BadVerbose', ...
        'verbose must be a logical scalar.');
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
        error('thole:solve_scf_sor:BadStopMetric', ...
            'stop_metric must be ''relres'' or ''max_dmu''.');
end

end

function tf = local_has_raw_rowcache_fast_path(op)

tf = isfield(op, 'row_cache') && ...
    isfield(op.row_cache, 'row_ptr') && ...
    isfield(op.row_cache, 'col_idx') && ...
    isfield(op.row_cache, 'dr') && ...
    isfield(op.row_cache, 'thole_f3') && ...
    isfield(op.row_cache, 'thole_f5') && ...
    isfield(op.row_cache, 'inv_r3_bare') && ...
    isfield(op.row_cache, 'inv_r5_bare');

if ~tf
    return;
end

rowPtr = op.row_cache.row_ptr(:);
colIdx = op.row_cache.col_idx(:);

if isempty(rowPtr)
    tf = false;
    return;
end

nEntries = rowPtr(end) - 1;

tf = numel(colIdx) == nEntries && ...
    size(op.row_cache.dr, 1) == nEntries && ...
    size(op.row_cache.dr, 2) == 3 && ...
    numel(op.row_cache.thole_f3) == nEntries && ...
    numel(op.row_cache.thole_f5) == nEntries && ...
    numel(op.row_cache.inv_r3_bare) == nEntries && ...
    numel(op.row_cache.inv_r5_bare) == nEntries;

end

function tf = local_should_compute_periodic_residual(iter, residualEvery)

if isempty(residualEvery) || residualEvery == 0 || isinf(residualEvery)
    tf = false;
    return;
end

% Match old-style diagnostics: compute on the first iteration and then at
% the requested cadence.
tf = (iter == 1) || (mod(iter, residualEvery) == 0);

end

function name = local_fast_path_name(tf)

if tf
    name = 'raw_rowcache';
else
    name = 'generic_apply_row';
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