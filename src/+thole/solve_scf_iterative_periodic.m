function [mu, scf] = solve_scf_iterative_periodic(sys, Eext, ewaldParams, scfParams)
%SOLVE_SCF_ITERATIVE_PERIODIC Matrix-free periodic fixed-point SCF solver.
%
% [mu, scf] = thole.solve_scf_iterative_periodic(sys, Eext, ewaldParams)
% [mu, scf] = thole.solve_scf_iterative_periodic(sys, Eext, ewaldParams, scfParams)
%
% This solver does NOT assemble the periodic dipole interaction matrix.
% Instead, each iteration applies the periodic operator through
% thole.induced_field_from_dipoles_thole_periodic(...).
%
% Inputs
%   sys         canonical polarization-system struct in atomic units
%   Eext        N x 3 external field in atomic units
%   ewaldParams struct with fields:
%                 .alpha
%                 .rcut
%                 .kcut
%                 .boundary   optional, default 'tinfoil'
%
%   scfParams   optional solver controls
%
% Supported controls in scfParams
%   .tol
%   .maxIter
%   .mixing
%   .verbose
%   .printEvery
%   .residualEvery
%   .stopMetric        'max_dmu' | 'relres'
%   .use_thole         logical, default true
%
% Output
%   mu         N x 3 induced dipoles
%   scf        convergence metadata
%
% Notes
%   - This is the periodic analog of the nonperiodic matrix-free iterative
%     solver, but it consumes periodic real-space and reciprocal-space caches.
%   - Relative residual is evaluated against the periodic fixed-point equation
%       mu = alpha .* (Eext + Edip(mu))
%     on the active polarizable sites only.

    if nargin < 4 || isempty(scfParams)
        scfParams = struct();
    end

    io.assert_atomic_units(sys);

    if nargin < 3 || isempty(ewaldParams)
        error('thole:solve_scf_iterative_periodic:MissingEwaldParams', ...
            'ewaldParams is required.');
    end
    if ~isfield(ewaldParams, 'alpha') || isempty(ewaldParams.alpha)
        error('thole:solve_scf_iterative_periodic:MissingAlpha', ...
            'ewaldParams.alpha is required.');
    end
    if ~isfield(ewaldParams, 'rcut') || isempty(ewaldParams.rcut)
        error('thole:solve_scf_iterative_periodic:MissingRcut', ...
            'ewaldParams.rcut is required.');
    end
    if ~isfield(ewaldParams, 'kcut') || isempty(ewaldParams.kcut)
        error('thole:solve_scf_iterative_periodic:MissingKcut', ...
            'ewaldParams.kcut is required.');
    end

    problem = thole.prepare_scf_problem(sys, Eext, scfParams);

    nSites        = problem.nSites;
    polMask       = problem.polMask;
    activeSites   = problem.activeSites(:);
    alpha         = problem.alpha;
    tol           = problem.tol;
    maxIter       = problem.maxIter;
    verbose       = problem.verbose;
    printEvery    = problem.printEvery;
    residualEvery = problem.residualEvery;

    mixing = 1.0;
    if isfield(problem, 'mixing') && ~isempty(problem.mixing)
        mixing = problem.mixing;
    elseif isfield(scfParams, 'mixing') && ~isempty(scfParams.mixing)
        mixing = scfParams.mixing;
    end

    stopMetric = "max_dmu";
    if isfield(problem, 'stopMetric') && ~isempty(problem.stopMetric)
        stopMetric = lower(string(problem.stopMetric));
    end

    use_thole = true;
    if isfield(scfParams, 'use_thole') && ~isempty(scfParams.use_thole)
        use_thole = logical(scfParams.use_thole);
    end

    % Build periodic caches once.
    cacheParams = struct();
    cacheParams.use_thole = use_thole;

    realCache = geom.build_periodic_realspace_cache(sys, problem, ewaldParams, cacheParams);
    kCache    = geom.build_periodic_kspace_cache(sys, problem, ewaldParams);

    mu = problem.mu0;

    history = zeros(maxIter, 1);
    relresHistory = nan(maxIter, 1);
    converged = false;

    dipoleParams = struct();
    dipoleParams.use_thole = use_thole;
    dipoleParams.problem = problem;
    dipoleParams.target_mask = polMask;
    dipoleParams.source_mask = polMask;
    dipoleParams.realspace_cache = realCache;
    dipoleParams.kspace_cache = kCache;

    if verbose
        fprintf(['SCF(iterative, periodic matrix-free): tol=%.3e, maxIter=%d, ' ...
                 'mixing=%.3f | nReal=%d | nK=%d\n'], ...
            tol, maxIter, mixing, realCache.nInteractions, kCache.num_kvec);
    end

    tSCF = tic;

    for iter = 1:maxIter
        Edip = thole.induced_field_from_dipoles_thole_periodic( ...
            sys, mu, ewaldParams, dipoleParams);

        Etot = Eext + Edip;

        mu_new = zeros(nSites, 3);
        mu_new(polMask, :) = alpha(polMask) .* Etot(polMask, :);

        mu_mixed = (1 - mixing) * mu + mixing * mu_new;
        delta = mu_mixed - mu;
        err = max(sqrt(sum(delta.^2, 2)));

        history(iter) = err;
        mu = mu_mixed;

        doResidual = strcmp(stopMetric, "relres") || ...
                     (iter == 1) || ...
                     (mod(iter, residualEvery) == 0);

        if doResidual
            relres = local_periodic_matrix_free_relres(sys, problem, mu, Eext, ewaldParams, dipoleParams);
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
            if ~strcmp(stopMetric, "relres") && isnan(relres)
                relres = local_periodic_matrix_free_relres(sys, problem, mu, Eext, ewaldParams, dipoleParams);
                relresHistory(iter) = relres;

                if verbose
                    fprintf(' final relres = %.3e\n', relres);
                end
            end

            history = history(1:iter);
            relresHistory = relresHistory(1:iter);

            if verbose
                fprintf('SCF(iterative, periodic matrix-free) converged in %d iterations (%.2f s)\n', ...
                    iter, toc(tSCF));
            end
            break;
        end
    end

    if ~converged
        history = history(1:maxIter);
        relresHistory = relresHistory(1:maxIter);

        if isnan(relresHistory(end))
            relresHistory(end) = local_periodic_matrix_free_relres(sys, problem, mu, Eext, ewaldParams, dipoleParams);
        end

        if verbose
            fprintf(['SCF(iterative, periodic matrix-free) hit maxIter=%d after %.2f s ' ...
                     '| final max|dmu| = %.3e\n'], ...
                maxIter, toc(tSCF), history(end));
        end
    end

    scf = struct();
    scf.method = 'iterative';
    scf.variant = 'periodic_matrix_free';
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
    scf.nRealInteractions = realCache.nInteractions;
    scf.num_kvec = kCache.num_kvec;
end

function relres = local_periodic_matrix_free_relres(sys, problem, mu, Eext, ewaldParams, dipoleParams)
%LOCAL_PERIODIC_MATRIX_FREE_RELRES Relative residual for periodic fixed-point SCF.
%
% Relative residual for:
%   mu = alpha .* (Eext + Edip(mu))
% on the active polarizable sites.

    polMask = problem.polMask;
    alpha = problem.alpha;

    Edip = thole.induced_field_from_dipoles_thole_periodic( ...
        sys, mu, ewaldParams, dipoleParams);

    rhs = zeros(problem.nSites, 3);
    rhs(polMask, :) = alpha(polMask) .* (Eext(polMask, :) + Edip(polMask, :));

    r = rhs(polMask, :) - mu(polMask, :);

    bn = norm(rhs(polMask, :), 'fro');
    if bn == 0
        bn = 1.0;
    end

    relres = norm(r, 'fro') / bn;
end