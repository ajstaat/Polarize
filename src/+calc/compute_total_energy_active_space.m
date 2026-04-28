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
% where all vectors are active-space stacked xyz vectors.
%
% The function also reports the stationary-form consistency check:
%
% At self-consistency,
%
%   A^{-1} mu = Eext + T mu
%
% therefore
%
%   E_pol = -0.5 * mu' * Eext
%
% In code:
%
%   external_charge_dipole = -mu' * Eext
%   total_stationary      = 0.5 * external_charge_dipole
%
% Inputs
%   sys
%       polarization system in atomic units
%
%   problem
%       struct from thole.prepare_scf_problem
%
%   mu
%       full-system N x 3 induced dipoles
%
%   Eext
%       full-system N x 3 external field
%
%   Tpol
%       dense active-space dipole interaction matrix
%
%   op
%       dense operator struct from thole.make_polarization_operator
%
% Output
%   energy struct with fields:
%       polarization_self
%       external_charge_dipole
%       dipole_dipole
%       total
%       total_stationary
%       stationary_consistency
%       relres
%       nPolSites
%       units

validate_inputs(sys, problem, mu, Eext);

Tpol = local_extract_Tpol_for_energy(opOrTpol);

io.assert_atomic_units(sys);

nVec = numel(problem.Eext_pol_vec);

if ~isequal(size(Tpol), [nVec nVec])
    error('calc:compute_total_energy_active_space:BadTpolSize', ...
        'Tpol must be %d x %d for this active-space problem.', nVec, nVec);
end

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

AinvVec = 1 ./ alphaVec;

polarizationSelf = 0.5 * sum((muVec.^2) .* AinvVec);
externalChargeDipole = -dot(muVec, EextVec);
dipoleDipole = -0.5 * dot(muVec, Tpol * muVec);

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

relres = thole.compute_active_space_relres(problem, Tpol, mu);

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

function Tpol = local_extract_Tpol_for_energy(opOrTpol)

if isnumeric(opOrTpol)
    Tpol = opOrTpol;
    return;
end

if ~isstruct(opOrTpol)
    error('calc:compute_total_energy_active_space:BadOperator', ...
        'Final input must be a dense Tpol matrix or an operator struct.');
end

if ~isfield(opOrTpol, 'kind') || ~strcmp(opOrTpol.kind, 'dense_matrix') || ...
        ~isfield(opOrTpol, 'Tpol') || isempty(opOrTpol.Tpol)

    backend = '<unknown>';
    if isfield(opOrTpol, 'backend') && ~isempty(opOrTpol.backend)
        backend = char(string(opOrTpol.backend));
    end

    error('calc:compute_total_energy_active_space:RequiresDenseOperator', ...
        ['Active-space energy currently requires a dense operator with op.Tpol ', ...
         'so it can evaluate mu''*T*mu. Requested backend "%s" is not dense. ', ...
         'Later matrix-free backends should provide an energy/apply-compatible path.'], backend);
end

Tpol = opOrTpol.Tpol;

if ~isnumeric(Tpol) || ndims(Tpol) ~= 2 || size(Tpol,1) ~= size(Tpol,2)
    error('calc:compute_total_energy_active_space:BadDenseOperatorMatrix', ...
        'op.Tpol must be a numeric square matrix.');
end

end