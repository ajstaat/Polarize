function [mu, scf] = solve_scf_iterative(sys, Eext, scfParams)
%SOLVE_SCF_ITERATIVE Solve induced dipoles by matrix-free fixed-point SCF iteration.
%
% This solver does NOT assemble the explicit dipole interaction matrix.
% Instead, each iteration applies the dipole-field operator through
% induced_field_from_dipoles_thole(...) or induced_field_from_dipoles(...).
%
% Intended role
%   - matrix-free iterative solver for larger systems
%   - useful when explicit dense Tpol assembly is too costly in memory
%
% Inputs
%   sys        canonical polarization-system struct in atomic units
%   Eext       N x 3 external field in atomic units
%   scfParams  optional solver controls
%
% Supported controls in scfParams
%   .tol
%   .maxIter
%   .mixing
%   .verbose
%   .printEvery
%   .residualEvery
%   .stopMetric        'max_dmu' | 'relres'
%   .softening
%   .use_thole
%
% Output
%   mu         N x 3 induced dipoles
%   scf        convergence metadata

    if nargin < 3 || isempty(scfParams)
        scfParams = struct();
    end

    io.assert_atomic_units(sys);

    problem = thole.prepare_scf_problem(sys, Eext, scfParams);

    nSites        = problem.nSites;
    polMask       = problem.polMask;
    alpha         = problem.alpha;

    tol           = problem.tol;
    maxIter       = problem.maxIter;
    mixing        = problem.mixing;
    verbose       = problem.verbose;
    printEvery    = problem.printEvery;
    residualEvery = problem.residualEvery;

    stopMetric = "max_dmu";
    if isfield(problem, 'stopMetric') && ~isempty(problem.stopMetric)
        stopMetric = lower(string(problem.stopMetric));
    end

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
    relresHistory = nan(maxIter, 1);
    converged = false;

    dipoleParams = struct();
    dipoleParams.exclude_self = true;
    dipoleParams.softening = softening;
    dipoleParams.target_mask = polMask;
    dipoleParams.source_mask = polMask;

    if isfield(scfParams, 'rcut') && ~isempty(scfParams.rcut)
        dipoleParams.rcut = scfParams.rcut;
    end

    if verbose
        fprintf('SCF(iterative, matrix-free): tol=%.3e, maxIter=%d, mixing=%.3f, use_thole=%d\n', ...
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

        doResidual = strcmp(stopMetric, "relres") || ...
                     (iter == 1) || (mod(iter, residualEvery) == 0) || (err < tol);

        if doResidual
            relres = local_matrix_free_relres(sys, problem, mu, Eext, use_thole, softening);
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

            if verbose
                fprintf('SCF(iterative, matrix-free) converged in %d iterations (%.2f s)\n', ...
                    iter, toc(tSCF));
            end
            break;
        end
    end

    if ~converged
        history = history(1:maxIter);
        relresHistory = relresHistory(1:maxIter);

        if verbose
            fprintf(['SCF(iterative, matrix-free) hit maxIter=%d after %.2f s ' ...
                     '| final max|dmu|=%.3e\n'], ...
                maxIter, toc(tSCF), history(end));
        end
    end

    scf = struct();
    scf.method = 'iterative';
    scf.variant = 'matrix_free';
    scf.converged = converged;
    scf.nIter = numel(history);
    scf.history = history;
    scf.relres_history = relresHistory;

    lastRelres = relresHistory(find(~isnan(relresHistory), 1, 'last'));
    if isempty(lastRelres)
        lastRelres = NaN;
    end
    scf.relres = lastRelres;
    scf.used_thole = use_thole;
    scf.used_matrix_solver = false;
end


function relres = local_matrix_free_relres(sys, problem, mu, Eext, use_thole, softening)
% Relative residual for the matrix-free fixed-point equation
%   mu = A * (Eext + Edip(mu))
% on the polarizable sites only.

    polMask = problem.polMask;
    alpha   = problem.alpha;

    dipoleParams = struct();
    dipoleParams.exclude_self = true;
    dipoleParams.softening = softening;
    dipoleParams.target_mask = polMask;
    dipoleParams.source_mask = polMask;

    if use_thole
        Edip = thole.induced_field_from_dipoles_thole(sys, mu, dipoleParams);
    else
        Edip = thole.induced_field_from_dipoles(sys, mu, dipoleParams);
    end

    rhs = zeros(problem.nSites, 3);
    rhs(polMask, :) = alpha(polMask) .* (Eext(polMask, :) + Edip(polMask, :));

    r = rhs(polMask, :) - mu(polMask, :);

    bn = norm(rhs(polMask, :), 'fro');
    if bn == 0
        bn = 1.0;
    end

    relres = norm(r, 'fro') / bn;
end