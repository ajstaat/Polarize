function energy = compute_total_energy_active_space(sys, problem, mu, Eext, opOrTpol)
%COMPUTE_TOTAL_ENERGY_ACTIVE_SPACE Compute active-space polarization energy.
%
% energy = calc.compute_total_energy_active_space(sys, problem, mu, Eext, Tpol)
% energy = calc.compute_total_energy_active_space(sys, problem, mu, Eext, op)
%
% Computes the full active-space polarization energy:
%
%   E_pol =
%       0.5 * mu' * A^{-1} * mu
%     -       mu' * Eext
%     - 0.5 * mu' * T * mu
%
% The final T*mu action may come from either:
%   - a dense active-space Tpol matrix
%   - an operator struct with op.apply(muVec)
%
% This allows large matrix-free workflows to compute energy without
% materializing dense Tpol.

validate_inputs(sys, problem, mu, Eext);

io.assert_atomic_units(sys);

activeSites = problem.activeSites(:);

muPol = mu(activeSites, :);
EextPol = Eext(activeSites, :);

muVec = util.stack_xyz(muPol);
EextVec = util.stack_xyz(EextPol);

alphaVec = problem.alpha_pol_vec(:);

if any(alphaVec <= 0)
    error('calc:compute_total_energy_active_space:NonpositiveAlpha', ...
        'All active-space alpha values must be positive for energy evaluation.');
end

nVec = numel(problem.Eext_pol_vec);

if numel(muVec) ~= nVec
    error('calc:compute_total_energy_active_space:BadMuVecSize', ...
        'Active-space mu vector length must equal problem.Eext_pol_vec length.');
end

TmuVec = local_apply_T(opOrTpol, muVec, nVec);

AinvVec = 1 ./ alphaVec;

polarizationSelf = 0.5 * sum((muVec.^2) .* AinvVec);
externalChargeDipole = -dot(muVec, EextVec);
dipoleDipole = -0.5 * dot(muVec, TmuVec);

total = polarizationSelf + externalChargeDipole + dipoleDipole;

% Stationary-form cross-check.
%
% At the SCF stationary point:
%
%   A^{-1} mu = Eext + T mu
%
% so:
%
%   E_pol =
%       0.5 mu' A^{-1} mu
%     -     mu' Eext
%     - 0.5 mu' T mu
%     = -0.5 mu' Eext
%
% Since externalChargeDipole = -mu' Eext:
%
%   E_stationary = 0.5 * externalChargeDipole
totalStationary = 0.5 * externalChargeDipole;
stationaryConsistency = total - totalStationary;

relres = local_compute_relres(problem, opOrTpol, mu, muVec, TmuVec);

energy = struct();
energy.polarization_self = polarizationSelf;
energy.external_charge_dipole = externalChargeDipole;
energy.dipole_dipole = dipoleDipole;
energy.total = total;

energy.total_stationary = totalStationary;
energy.stationary_consistency = stationaryConsistency;

energy.relres = relres;
energy.nPolSites = problem.nPolSites;

energy.units = struct();
energy.units.energy = 'hartree';
energy.units.length = 'bohr';
energy.units.alpha = 'atomic_unit';
energy.units.charge = 'elementary_charge';

end

% =========================================================================
% Operator helpers
% =========================================================================

function TmuVec = local_apply_T(opOrTpol, muVec, nVec)

if isnumeric(opOrTpol)
    Tpol = opOrTpol;

    if ~isequal(size(Tpol), [nVec nVec])
        error('calc:compute_total_energy_active_space:BadTpolSize', ...
            'Tpol must be %d x %d for this active-space problem.', nVec, nVec);
    end

    TmuVec = Tpol * muVec;
    return;
end

if ~isstruct(opOrTpol)
    error('calc:compute_total_energy_active_space:BadOperator', ...
        'Final input must be a dense Tpol matrix or an operator struct.');
end

if isfield(opOrTpol, 'kind') && strcmp(opOrTpol.kind, 'dense_matrix')
    if ~isfield(opOrTpol, 'Tpol') || isempty(opOrTpol.Tpol)
        error('calc:compute_total_energy_active_space:MissingDenseMatrix', ...
            'Dense operator must contain op.Tpol.');
    end

    Tpol = opOrTpol.Tpol;

    if ~isequal(size(Tpol), [nVec nVec])
        error('calc:compute_total_energy_active_space:BadDenseOperatorMatrixSize', ...
            'op.Tpol must be %d x %d for this active-space problem.', nVec, nVec);
    end

    TmuVec = Tpol * muVec;
    return;
end

if isfield(opOrTpol, 'apply') && isa(opOrTpol.apply, 'function_handle')
    TmuVec = opOrTpol.apply(muVec);

    if numel(TmuVec) ~= nVec
        error('calc:compute_total_energy_active_space:BadApplyOutputSize', ...
            'op.apply(muVec) must return a vector of length %d.', nVec);
    end

    TmuVec = TmuVec(:);
    return;
end

backend = '<unknown>';
if isfield(opOrTpol, 'backend') && ~isempty(opOrTpol.backend)
    backend = char(string(opOrTpol.backend));
end

error('calc:compute_total_energy_active_space:OperatorMissingApply', ...
    ['Energy evaluation requires either dense op.Tpol or op.apply(muVec). ', ...
     'Requested backend "%s" provides neither.'], backend);

end

function relres = local_compute_relres(problem, opOrTpol, mu, muVec, TmuVec)

alphaVec = problem.alpha_pol_vec(:);
Evec = problem.Eext_pol_vec(:);

resVec = muVec - alphaVec .* (Evec + TmuVec);

rhsScale = norm(alphaVec .* Evec);
if rhsScale == 0
    rhsScale = 1;
end

relres = norm(resVec) / rhsScale;

% For dense/raw Tpol, this should match the canonical routine. Do not call
% thole.compute_active_space_relres here because that routine currently
% expects dense Tpol and would defeat matrix-free workflows.
if nargin < 5 %#ok<UNRCH>
    relres = thole.compute_active_space_relres(problem, opOrTpol, mu);
end

end

% =========================================================================
% Validation helpers
% =========================================================================

function validate_inputs(sys, problem, mu, Eext)

if ~isstruct(sys)
    error('calc:compute_total_energy_active_space:BadSys', ...
        'sys must be a struct.');
end

if ~isstruct(problem)
    error('calc:compute_total_energy_active_space:BadProblem', ...
        'problem must be a struct from thole.prepare_scf_problem.');
end

requiredSys = {
    'site_pos'
    'site_alpha'
    'site_is_polarizable'
    'units'
};

for k = 1:numel(requiredSys)
    name = requiredSys{k};

    if ~isfield(sys, name) || isempty(sys.(name))
        error('calc:compute_total_energy_active_space:MissingSysField', ...
            'sys.%s is required and missing/empty.', name);
    end
end

requiredProblem = {
    'activeSites'
    'nPolSites'
    'Eext_pol_vec'
    'alpha_pol_vec'
};

for k = 1:numel(requiredProblem)
    name = requiredProblem{k};

    if ~isfield(problem, name) || isempty(problem.(name))
        error('calc:compute_total_energy_active_space:MissingProblemField', ...
            'problem.%s is required and missing/empty.', name);
    end
end

nSites = size(sys.site_pos, 1);

if ~isnumeric(mu) || ~isequal(size(mu), [nSites 3])
    error('calc:compute_total_energy_active_space:BadMu', ...
        'mu must be nSites x 3.');
end

if ~isnumeric(Eext) || ~isequal(size(Eext), [nSites 3])
    error('calc:compute_total_energy_active_space:BadEext', ...
        'Eext must be nSites x 3.');
end

if numel(problem.activeSites) ~= problem.nPolSites
    error('calc:compute_total_energy_active_space:BadActiveSites', ...
        'numel(problem.activeSites) must equal problem.nPolSites.');
end

nVec = 3 * problem.nPolSites;

if numel(problem.Eext_pol_vec) ~= nVec
    error('calc:compute_total_energy_active_space:BadEextPolVec', ...
        'problem.Eext_pol_vec must have length 3*problem.nPolSites.');
end

if numel(problem.alpha_pol_vec) ~= nVec
    error('calc:compute_total_energy_active_space:BadAlphaPolVec', ...
        'problem.alpha_pol_vec must have length 3*problem.nPolSites.');
end

end