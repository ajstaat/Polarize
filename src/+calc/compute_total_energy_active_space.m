function energy = compute_total_energy_active_space(sys, problem, mu, Eext, Tpol)
%COMPUTE_TOTAL_ENERGY_ACTIVE_SPACE Total polarization energy using Tpol.
%
% Returns the same fields as the legacy total-energy routine.
%
% Notes
%   - Computation is performed in the polarizable-site active space.
%   - The reported total uses the stationary induced-dipole expression:
%         E = 0.5 * (-mu·Eext) + (-0.5 * mu·Tmu)

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

    mu_pol       = mu(activeSites, :);
    mu_pol_vec   = util.stack_xyz(mu_pol);
    Eext_pol     = Eext(activeSites, :);
    Eext_pol_vec = util.stack_xyz(Eext_pol);

    % External-field contribution
    e_ext = -dot(mu_pol_vec, Eext_pol_vec);

    % Dipole-dipole contribution
    e_dd = -0.5 * dot(mu_pol_vec, Tpol * mu_pol_vec);

    energy = struct();
    energy.polarization_self = NaN;
    energy.external_charge_dipole = e_ext;
    energy.dipole_dipole = e_dd;
    energy.total = 0.5 * e_ext + e_dd;
end