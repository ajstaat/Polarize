function energy = compute_total_energy_active_space(sys, problem, mu, Eext, Tpol)
%COMPUTE_TOTAL_ENERGY_ACTIVE_SPACE Total polarization energy using Tpol.
%
% Returns the same fields as the legacy total-energy routine.
%
% Notes
%   - Computation is performed in the polarizable-site active space.
%   - The reported decomposition is:
%         polarization_self      = 0.5 * mu' * Ainv * mu
%         external_charge_dipole = -mu' * Eext
%         dipole_dipole          = -0.5 * mu' * Tpol * mu
%
%     so that
%         total = polarization_self + external_charge_dipole + dipole_dipole
%
%   - At the stationary point, this is equivalent to
%         total = 0.5 * (-mu·Eext) + (-0.5 * mu·Tmu)

    %#ok<INUSD> sys is kept for API compatibility with legacy/full-space routines

    activeSites = problem.activeSites(:);
    nPolSites   = problem.nPolSites;

    if size(mu, 2) ~= 3
        error('mu must be N x 3.');
    end
    if size(Eext, 2) ~= 3
        error('Eext must be N x 3.');
    end
    if numel(activeSites) ~= nPolSites
        error('problem.activeSites and problem.nPolSites are inconsistent.');
    end
    if ~isequal(size(Tpol), [3*nPolSites, 3*nPolSites])
        error('Tpol must be 3*Np x 3*Np for the active polarizable subspace.');
    end
    if ~isfield(problem, 'alpha_pol_vec') || isempty(problem.alpha_pol_vec)
        error('compute_total_energy_active_space requires problem.alpha_pol_vec.');
    end

    mu_pol       = mu(activeSites, :);
    mu_pol_vec   = util.stack_xyz(mu_pol);
    Eext_pol     = Eext(activeSites, :);
    Eext_pol_vec = util.stack_xyz(Eext_pol);

    alpha_pol_vec = problem.alpha_pol_vec(:);
    if numel(alpha_pol_vec) ~= numel(mu_pol_vec)
        error('problem.alpha_pol_vec has inconsistent length.');
    end
    if any(alpha_pol_vec == 0)
        error('problem.alpha_pol_vec contains zero entries.');
    end

    % Self / inverse-polarizability contribution
    e_self = 0.5 * dot(mu_pol_vec, mu_pol_vec ./ alpha_pol_vec);

    % External-field contribution
    e_ext = -dot(mu_pol_vec, Eext_pol_vec);

    % Dipole-dipole contribution
    e_dd = -0.5 * dot(mu_pol_vec, Tpol * mu_pol_vec);

    energy = struct();
    energy.polarization_self = e_self;
    energy.external_charge_dipole = e_ext;
    energy.dipole_dipole = e_dd;
    energy.total = e_self + e_ext + e_dd;

    % Optional stationary-form cross-check
    energy.total_stationary = 0.5 * e_ext + e_dd;
    energy.stationary_consistency = energy.total - energy.total_stationary;
end