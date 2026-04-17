function [mu, scf] = solve_scf_sor(problem, Tpol)
%SOLVE_SCF_SOR Solve induced dipoles by block GS / SOR sweeps on Tpol.
%
% Uses the fixed-point form on the polarizable-only active space:
%   mu_i = alpha_i * (Eext_i + sum_j T_ij mu_j)
%
% with Gauss-Seidel / SOR ordering:
%   earlier j use newest values
%   later   j use previous-iteration values
%
% omega = 1.0 -> GS
% 0 < omega < 1 -> under-relaxed GS
% 1 < omega < 2 -> SOR

    nSites    = problem.nSites;
    nPolSites = problem.nPolSites;

    if ~isequal(size(Tpol), [3*nPolSites, 3*nPolSites])
        error('Tpol must be 3*Np x 3*Np for the active polarizable subspace.');
    end

    activeSites   = problem.activeSites(:);
    alpha_pol     = problem.alpha_pol(:);
    Eext_pol      = problem.Eext_pol;
    mu_pol        = problem.mu0_pol;

    if numel(activeSites) ~= nPolSites
        error('problem.activeSites and problem.nPolSites are inconsistent.');
    end
    if numel(alpha_pol) ~= nPolSites
        error('problem.alpha_pol must have length nPolSites.');
    end
    if size(Eext_pol, 1) ~= nPolSites || size(Eext_pol, 2) ~= 3
        error('problem.Eext_pol must be nPolSites x 3.');
    end
    if size(mu_pol, 1) ~= nPolSites || size(mu_pol, 2) ~= 3
        error('problem.mu0_pol must be nPolSites x 3.');
    end

    tol           = problem.tol;
    maxIter       = problem.maxIter;
    omega         = problem.omega;
    verbose       = problem.verbose;
    printEvery    = problem.printEvery;
    residualEvery = problem.residualEvery;

    stopMetric = "max_dmu";
    if isfield(problem, 'stopMetric') && ~isempty(problem.stopMetric)
        stopMetric = lower(string(problem.stopMetric));
    end

    if omega <= 0 || omega >= 2
        error('problem.omega must satisfy 0 < omega < 2.');
    end

    bpol = problem.alpha_pol_vec(:) .* problem.Eext_pol_vec(:);
    bn = norm(bpol);
    if bn == 0
        bn = 1.0;
    end

    history = zeros(maxIter,1);
    relresHistory = nan(maxIter,1);
    converged = false;

    if verbose
        if abs(omega - 1.0) < 1e-14
            solverName = 'GS';
        else
            solverName = 'SOR';
        end
        fprintf('SCF(%s): tol=%.3e, maxIter=%d, omega=%.3f, nPol=%d\n', ...
            solverName, tol, maxIter, omega, nPolSites);
    end

    for iter = 1:maxIter
        mu_old_pol = mu_pol;

        % Stack once per sweep so we can use block matvecs
        mu_pol_vec     = util.stack_xyz(mu_pol);
        mu_old_pol_vec = util.stack_xyz(mu_old_pol);

        for k = 1:nPolSites
            Ik = (3*k-2):(3*k);

            E_loc = Eext_pol(k,:).';

            % Earlier active sites: newest values
            if k > 1
                Jprev = 1:(3*(k-1));
                E_loc = E_loc + Tpol(Ik, Jprev) * mu_pol_vec(Jprev);
            end

            % Later active sites: previous-iteration values
            if k < nPolSites
                Jnext = (3*k+1):(3*nPolSites);
                E_loc = E_loc + Tpol(Ik, Jnext) * mu_old_pol_vec(Jnext);
            end

            mu_gs_k = alpha_pol(k) * E_loc;
            mu_pol(k,:) = ((1 - omega) * mu_old_pol(k,:).' + omega * mu_gs_k).';

            % Keep stacked "newest values" vector in sync for subsequent rows
            mu_pol_vec(Ik) = mu_pol(k,:).';
        end

        dmu_pol = mu_pol - mu_old_pol;
        err = max(sqrt(sum(dmu_pol.^2, 2)));
        history(iter) = err;

        doResidual = strcmp(stopMetric, "relres") || ...
                     (iter == 1) || (mod(iter, residualEvery) == 0) || (err < tol);

        if doResidual
            relres = thole.compute_active_space_relres(problem, Tpol, mu_pol);
            relresHistory(iter) = relres;
        else
            relres = NaN;
        end

        if verbose && (iter == 1 || mod(iter, printEvery) == 0 || err < tol)
            if isnan(relres)
                fprintf(' iter %4d | max|dmu| = %.3e\n', iter, err);
            else
                fprintf(' iter %4d | max|dmu| = %.3e | relres = %.3e\n', iter, err, relres);
            end
        end

        switch stopMetric
            case "relres"
                if ~isnan(relres) && relres < tol
                    converged = true;
                end
            otherwise
                if err < tol
                    converged = true;
                end
        end

        if converged
            history = history(1:iter);
            relresHistory = relresHistory(1:iter);
            break;
        end
    end

    if ~converged
        history = history(1:maxIter);
        relresHistory = relresHistory(1:maxIter);
    end

    mu = zeros(nSites, 3);
    if nPolSites > 0
        mu(activeSites, :) = mu_pol;
    end

    scf = struct();
    if abs(omega - 1.0) < 1e-14
        scf.method = 'gs';
    else
        scf.method = 'sor';
    end
    scf.omega = omega;
    scf.converged = converged;
    scf.nIter = numel(history);
    scf.history = history;
    scf.relres_history = relresHistory;

    lastRelres = relresHistory(find(~isnan(relresHistory), 1, 'last'));
    if isempty(lastRelres)
        lastRelres = NaN;
    end
    scf.relres = lastRelres;
    scf.used_matrix_solver = true;
    scf.nPolSites = nPolSites;
end