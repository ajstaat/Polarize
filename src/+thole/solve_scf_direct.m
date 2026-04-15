function [mu, direct] = solve_scf_direct(problem, Tpol)
%SOLVE_SCF_DIRECT Solve induced dipoles by direct linear solve in active space.
%
% Solves
%   (I - A*T) mu = A*Eext
%
% on the polarizable-only active space, then expands the result back to the
% full N x 3 dipole array.

    nSites    = problem.nSites;
    nPolSites = problem.nPolSites;

    if ~isequal(size(Tpol), [3*nPolSites, 3*nPolSites])
        error('Tpol must be 3*Np x 3*Np for the active polarizable subspace.');
    end

    activeSites    = problem.activeSites;
    alpha_pol_vec  = problem.alpha_pol_vec;
    Eext_pol_vec   = problem.Eext_pol_vec;

    % Build active-space linear system:
    %   (I - A*T) mu = A*Eext
    Ipol = eye(3*nPolSites);
    Mpol = Ipol - (alpha_pol_vec .* Tpol);
    bpol = alpha_pol_vec .* Eext_pol_vec;

    mu_pol_vec = Mpol \ bpol;

    % Expand back to full N x 3 array
    mu = zeros(nSites, 3);
    if nPolSites > 0
        mu(activeSites, :) = util.unstack_xyz(mu_pol_vec);
    end

    % Diagnostics
    r = bpol - Mpol * mu_pol_vec;
    residual_norm = norm(r);
    relres = thole.compute_active_space_relres(problem, Tpol, mu);

    direct = struct();
    direct.method = 'direct';
    direct.used_matrix_solver = true;
    direct.nPolSites = nPolSites;
    direct.residual_norm = residual_norm;
    direct.relres = relres;
end