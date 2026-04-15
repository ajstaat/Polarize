function [mu, scf] = solve_scf_matrix_iterative(sys, Eext, T, scfParams)
%SOLVE_SCF_MATRIX_ITERATIVE Solve SCF using explicit interaction matrix T.
%
% Fixed-point / Jacobi-style update with optional linear mixing:
%   mu <- (1-mixing)*mu + mixing*A*(Eext + T*mu)

    problem = thole.prepare_scf_problem(sys, Eext, scfParams);

    nSites = problem.nSites;
    if ~isequal(size(T), [3*nSites, 3*nSites])
        error('T must be 3N x 3N.');
    end

    polMask = problem.polMask;
    alpha   = problem.alpha;

    mu_vec  = problem.mu0_vec;
    Eext_vec = problem.Eext_vec;

    tol       = problem.tol;
    maxIter   = problem.maxIter;
    mixing    = problem.mixing;
    verbose   = problem.verbose;
    printEvery = problem.printEvery;

    history = zeros(maxIter, 1);
    converged = false;

    if verbose
        fprintf('SCF(matrix_iterative): tol=%.3e, maxIter=%d, mixing=%.3f\n', ...
            tol, maxIter, mixing);
    end

    tSCF = tic;

    for iter = 1:maxIter
        Edip_vec = T * mu_vec;
        Etot_vec = Eext_vec + Edip_vec;
        Etot = util.unstack_xyz(Etot_vec);

        mu_new = zeros(nSites, 3);
        mu_new(polMask, :) = alpha(polMask) .* Etot(polMask, :);

        mu_new_vec = util.stack_xyz(mu_new);
        mu_mixed_vec = (1 - mixing) * mu_vec + mixing * mu_new_vec;

        delta_vec = mu_mixed_vec - mu_vec;
        delta = util.unstack_xyz(delta_vec);
        err = max(sqrt(sum(delta.^2, 2)));

        history(iter) = err;
        mu_vec = mu_mixed_vec;

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

    mu = util.unstack_xyz(mu_vec);

    scf = struct();
    scf.method = 'matrix_iterative';
    scf.converged = converged;
    scf.nIter = numel(history);
    scf.history = history;
    scf.used_matrix_solver = true;
end