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

    mu_pol_vec = util.stack_xyz(mu_pol);

    alpha_pol = problem.alpha_pol;
    Eext_pol  = problem.Eext_pol;

    bpol = zeros(3*nPolSites, 1);
    for k = 1:nPolSites
        Ik = (3*k-2):(3*k);
        bpol(Ik) = alpha_pol(k) * Eext_pol(k,:).';
    end

    Tmu_pol = Tpol * mu_pol_vec;
    r = bpol - mu_pol_vec;

    for k = 1:nPolSites
        Ik = (3*k-2):(3*k);
        r(Ik) = r(Ik) + alpha_pol(k) * Tmu_pol(Ik);
    end

    bn = norm(bpol);
    if bn == 0
        bn = 1.0;
    end

    relres = norm(r) / bn;
end