function relres = compute_active_space_relres(problem, Tpol, mu)
%COMPUTE_ACTIVE_SPACE_RELRES Relative residual for (I - A*T) mu = A*Eext
% on the polarizable-only active space.
%
% Input mu may be either:
%   - full-space [N x 3]
%   - active-space [Np x 3]

    nPolSites = problem.nPolSites;

    if isequal(size(mu), [problem.nSites, 3])
        mu_pol = mu(problem.activeSites, :);
    elseif isequal(size(mu), [nPolSites, 3])
        mu_pol = mu;
    else
        error('mu must be either N x 3 or Np x 3.');
    end

    if ~isequal(size(Tpol), [3*nPolSites, 3*nPolSites])
        error('Tpol must be 3*Np x 3*Np.');
    end

    mu_pol_vec    = util.stack_xyz(mu_pol);
    alpha_pol_vec = problem.alpha_pol_vec(:);
    Eext_pol_vec  = problem.Eext_pol_vec(:);

    if numel(alpha_pol_vec) ~= 3*nPolSites || numel(Eext_pol_vec) ~= 3*nPolSites
        error('problem alpha/Eext active-space vectors are inconsistent with nPolSites.');
    end

    bpol = alpha_pol_vec .* Eext_pol_vec;
    r = bpol - mu_pol_vec + alpha_pol_vec .* (Tpol * mu_pol_vec);

    bn = norm(bpol);
    if bn == 0
        bn = 1.0;
    end

    relres = norm(r) / bn;
end