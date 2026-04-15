function [mu, scf] = solve_scf_iterative(sys, Eext, scfParams)
%SOLVE_SCF_ITERATIVE Solve induced dipoles by fixed-point SCF iteration.
%
% Legacy pairwise-field implementation. Kept mainly for validation/debug.
% For performance, prefer solve_scf_matrix_iterative or solve_scf_sor.

    problem = thole.prepare_scf_problem(sys, Eext, scfParams);

    nSites = problem.nSites;
    polMask = problem.polMask;
    alpha   = problem.alpha;

    tol       = problem.tol;
    maxIter   = problem.maxIter;
    mixing    = problem.mixing;
    verbose   = problem.verbose;
    printEvery = problem.printEvery;

    mu = problem.mu0;

    softening = 0.0;
    if isfield(scfParams, 'softening') && ~isempty(scfParams.softening)
        softening = scfParams.softening;
    end

    use_thole = true;
    if isfield(scfParams, 'use_thole') && ~isempty(scfParams.use_thole)
        use_thole = logical(scfParams.use_thole);
    end

    history = zeros(maxIter, 1);
    converged = false;

    dipoleParams = struct();
    dipoleParams.exclude_self = true;
    dipoleParams.softening = softening;
    dipoleParams.target_mask = polMask;
    dipoleParams.source_mask = polMask;

    if verbose
        fprintf('SCF(iterative): tol=%.3e, maxIter=%d, mixing=%.3f, use_thole=%d\n', ...
            tol, maxIter, mixing, use_thole);
    end

    tSCF = tic;

    for iter = 1:maxIter
        if use_thole
            Edip = thole.induced_field_from_dipoles_thole(sys, mu, dipoleParams);
        else
            Edip = thole.induced_field_from_dipoles(sys, mu, dipoleParams);
        end

        Etot = Eext + Edip;

        mu_new = zeros(nSites, 3);
        mu_new(polMask, :) = alpha(polMask) .* Etot(polMask, :);

        mu_mixed = (1 - mixing) * mu + mixing * mu_new;
        delta = mu_mixed - mu;
        err = max(sqrt(sum(delta.^2, 2)));

        history(iter) = err;
        mu = mu_mixed;

        if verbose && (iter == 1 || mod(iter, printEvery) == 0 || err < tol)
            fprintf(' iter %4d | max|dmu| = %.3e\n', iter, err);
        end

        if err < tol
            converged = true;
            history = history(1:iter);

            if verbose
                fprintf('SCF(iterative) converged in %d iterations (%.2f s)\n', ...
                    iter, toc(tSCF));
            end
            break;
        end
    end

    if ~converged
        history = history(1:maxIter);
        if verbose
            fprintf('SCF(iterative) hit maxIter=%d after %.2f s | final max|dmu|=%.3e\n', ...
                maxIter, toc(tSCF), history(end));
        end
    end

    scf = struct();
    scf.method = 'iterative';
    scf.converged = converged;
    scf.nIter = numel(history);
    scf.history = history;
    scf.used_thole = use_thole;
end