function energy = compute_total_energy_active_space(sys, problem, mu, Eext, Tpol)
%COMPUTE_TOTAL_ENERGY_ACTIVE_SPACE Total polarization energy using Tpol.
%
% Returns the same fields as the legacy total-energy routine.

    mu_pol = mu(problem.activeSites, :);
    mu_pol_vec = util.stack_xyz(mu_pol);
    Eext_pol   = Eext(problem.activeSites, :);
    Eext_pol_vec = util.stack_xyz(Eext_pol);

    % Self + external-field contribution
    e_ext = -dot(mu_pol_vec, Eext_pol_vec);

    % Dipole-dipole contribution
    e_dd = -0.5 * dot(mu_pol_vec, Tpol * mu_pol_vec);

    % Keep the existing naming convention as closely as possible
    energy = struct();
    energy.polarization_self = NaN;   % split not uniquely available unless you mirror legacy convention exactly
    energy.external_charge_dipole = e_ext;
    energy.dipole_dipole = e_dd;
    energy.total = 0.5 * e_ext + e_dd;
end