function [mu, scf] = solve_scf_matrix_iterative(problem, Tpol)
%SOLVE_SCF_MATRIX_ITERATIVE Solve SCF using explicit active-space matrix Tpol.
%
% Fixed-point / Jacobi-style update with optional linear mixing:
%   mu_pol <- (1-mixing)*mu_pol + mixing*(alpha_pol .* (Eext_pol + Tpol*mu_pol))
%
% Inputs
%   problem  struct from thole.prepare_scf_problem(...)
%   Tpol     [3*Np x 3*Np] polarizable-only interaction matrix
%
% Output
%   mu       [N x 3] full-space induced dipoles (zeros on nonpolarizable sites)
%   scf      struct with convergence metadata

    nSites   = problem.nSites;
    nPolSites = problem.nPolSites;

    if ~isequal(size(Tpol), [3*nPolSites, 3*nPolSites])
        error('Tpol must be 3*Np x 3*Np for the active polarizable subspace.');
    end

    activeSites    = problem.activeSites;
    alpha_pol_vec  = problem.alpha_pol_vec;
    Eext_pol_vec   = problem.Eext_pol_vec;
    mu_pol_vec     = problem.mu0_pol_vec;

    tol          = problem.tol;
    maxIter      = problem.maxIter;
    mixing       = problem.mixing;
    verbose      = problem.verbose;
    printEvery   = problem.printEvery;
    stopMetric   = 'max_dmu';

    if isfield(problem, 'stopMetric') && ~isempty(problem.stopMetric)
        stopMetric = lower(string(problem.stopMetric));
    end

    history = zeros(maxIter, 1);
    converged = false;

    if verbose
        fprintf('SCF(matrix_iterative): tol=%.3e, maxIter=%d, mixing=%.3f, nPol=%d\n', ...
            tol, maxIter, mixing, nPolSites);
    end

    tSCF = tic;

    for iter = 1:maxIter
        % induced field from other dipoles, in active space
        Edip_pol_vec = Tpol * mu_pol_vec;

        % total field on polarizable sites
        Etot_pol_vec = Eext_pol_vec + Edip_pol_vec;

        % Jacobi/fixed-point update: mu = alpha .* E
        mu_new_pol_vec = alpha_pol_vec .* Etot_pol_vec;

        % linear mixing
        mu_mixed_pol_vec = (1 - mixing) * mu_pol_vec + mixing * mu_new_pol_vec;

        % convergence metric
        delta_pol_vec = mu_mixed_pol_vec - mu_pol_vec;

        switch stopMetric
            case "max_dmu"
                delta_pol = util.unstack_xyz(delta_pol_vec);
                err = max(sqrt(sum(delta_pol.^2, 2)));

            otherwise
                delta_pol = util.unstack_xyz(delta_pol_vec);
                err = max(sqrt(sum(delta_pol.^2, 2)));
        end

        history(iter) = err;
        mu_pol_vec = mu_mixed_pol_vec;

        if verbose && (iter == 1 || mod(iter, printEvery) == 0 || err < tol)
            fprintf(' iter %4d | max|dmu| = %.3e\n', iter, err);
        end

        if err < tol
            converged = true;
            history = history(1:iter);

            if verbose
                fprintf('SCF(matrix_iterative) converged in %d iterations (%.2f s)\n', ...
                    iter, toc(tSCF));
            end
            break;
        end
    end

    if ~converged
        history = history(1:maxIter);
        if verbose
            fprintf('SCF(matrix_iterative) hit maxIter=%d after %.2f s | final max|dmu|=%.3e\n', ...
                maxIter, toc(tSCF), history(end));
        end
    end

    % Expand active-space solution back to full N x 3 output
    mu = zeros(nSites, 3);
    if nPolSites > 0
        mu(activeSites, :) = util.unstack_xyz(mu_pol_vec);
    end

    scf = struct();
    scf.method = 'matrix_iterative';
    scf.converged = converged;
    scf.nIter = numel(history);
    scf.history = history;
    scf.used_matrix_solver = true;
    scf.nPolSites = nPolSites;
end