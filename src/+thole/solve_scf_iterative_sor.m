function [mu, scf] = solve_scf_iterative_sor(sys, Eext, scfParams)
%SOLVE_SCF_ITERATIVE_SOR Matrix-free GS/SOR solve using row-wise cached neighbors.
%
% This solver does NOT assemble dense Tpol. It performs row sweeps over a
% directed active-space neighbor cache derived from the nonperiodic pair cache.
%
% Fixed-point update on active polarizable sites:
%   mu_i = alpha_i * (Eext_i + sum_j T_ij mu_j)
%
% SOR update:
%   mu_i <- (1-omega) * mu_i_old + omega * mu_i_GS
%
% Inputs
%   sys        canonical polarization-system struct in atomic units
%   Eext       N x 3 external field
%   scfParams  optional solver controls
%
% Supported controls in scfParams
%   .tol
%   .maxIter
%   .omega
%   .verbose
%   .printEvery
%   .residualEvery
%   .stopMetric        'max_dmu' | 'relres'
%   .softening
%   .rcut
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
    activeSites   = problem.activeSites(:);
    nPolSites     = problem.nPolSites;
    alpha_pol     = problem.alpha_pol(:);
    Eext_pol      = problem.Eext_pol;

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
        error('thole:solve_scf_iterative_sor:BadOmega', ...
            'scfParams.omega must satisfy 0 < omega < 2.');
    end

    softening = 0.0;
    if isfield(scfParams, 'softening') && ~isempty(scfParams.softening)
        softening = scfParams.softening;
    end

    % Build unordered pair cache once, then directed row cache.
    cacheOpts = struct();
    cacheOpts.site_mask = polMask;

    if isfield(scfParams, 'rcut') && ~isempty(scfParams.rcut)
        cacheOpts.rcut = scfParams.rcut;
    else
        cacheOpts.rcut = inf;
    end

    pairCache = geom.build_nonperiodic_pair_cache(sys, cacheOpts);
    rowCache = geom.build_active_row_cache(sys, problem, pairCache);

    mu = problem.mu0;
    mu_pol = mu(activeSites, :);

    history = zeros(maxIter, 1);
    relresHistory = nan(maxIter, 1);
    converged = false;

    % For residual checks, reuse the fast symmetric full-apply cache.
    dipoleParams = struct();
    dipoleParams.exclude_self = true;
    dipoleParams.softening = softening;
    dipoleParams.target_mask = polMask;
    dipoleParams.source_mask = polMask;
    if isfield(scfParams, 'rcut') && ~isempty(scfParams.rcut)
        dipoleParams.rcut = scfParams.rcut;
    end
    dipoleParams.geom_cache = pairCache;

    % Pull row-cache arrays locally once.
    row_ptr = rowCache.row_ptr;
    col_idx = rowCache.col_idx;
    dr_all  = rowCache.dr;

    haveThole = isfield(rowCache, 'thole_f3') && isfield(rowCache, 'thole_f5');
    if haveThole
        f3_all = rowCache.thole_f3;
        f5_all = rowCache.thole_f5;
    else
        % This fallback should almost never be needed if pair cache was built
        % from a normal polarization system.
        r_bare_all = rowCache.r_bare;
        alpha_full = sys.site_alpha(:);
        a = sys.thole_a;
    end

    if softening == 0
        useSoftening = false;
        invR3_all = rowCache.inv_r3_bare;
        invR5_all = rowCache.inv_r5_bare;
        r2_bare_all = [];
    else
        useSoftening = true;
        invR3_all = [];
        invR5_all = [];
        r2_bare_all = rowCache.r2_bare;
    end

    if verbose
        fprintf('SCF(iterative, matrix-free SOR): tol=%.3e, maxIter=%d, omega=%.3f\n', ...
            tol, maxIter, omega);
    end

    tSCF = tic;

    for iter = 1:maxIter
        mu_old_pol = mu_pol;

        for i = 1:nPolSites
            k0 = row_ptr(i);
            k1 = row_ptr(i + 1) - 1;

            E_loc = Eext_pol(i, :);

            if k1 >= k0
                idx = k0:k1;
                cols = col_idx(idx);

                muNbr = mu_pol(cols, :);      % current in-place values: GS/SOR ordering handled naturally
                dr    = dr_all(idx, :);       % directed r_j - r_i for target i <- source j

                if useSoftening
                    r2 = r2_bare_all(idx) + softening^2;
                    invR = 1 ./ sqrt(r2);
                    invR3 = invR ./ r2;
                    invR5 = invR3 ./ r2;
                else
                    invR3 = invR3_all(idx);
                    invR5 = invR5_all(idx);
                end

                if haveThole
                    f3 = f3_all(idx);
                    f5 = f5_all(idx);
                else
                    % Rare fallback path
                    nNbr = numel(idx);
                    f3 = zeros(nNbr, 1);
                    f5 = zeros(nNbr, 1);
                    iFull = activeSites(i);
                    jFull = activeSites(cols);
                    for kk = 1:nNbr
                        tf = thole.thole_f3f5_factors(r_bare_all(idx(kk)), ...
                            alpha_full(iFull), alpha_full(jFull(kk)), a);
                        f3(kk) = tf.f3;
                        f5(kk) = tf.f5;
                    end
                end

                muDotR = sum(muNbr .* dr, 2);

                coeff1 = 3 .* (f5 .* muDotR .* invR5);   % nNbr x 1
                coeff2 = (f3 .* invR3);                  % nNbr x 1

                contrib = coeff1 .* dr - coeff2 .* muNbr;   % nNbr x 3
                E_loc = E_loc + sum(contrib, 1);
            end

            mu_gs_i = alpha_pol(i) .* E_loc;
            mu_pol(i, :) = (1 - omega) .* mu_old_pol(i, :) + omega .* mu_gs_i;
        end

        dmu_pol = mu_pol - mu_old_pol;
        err = max(sqrt(sum(dmu_pol.^2, 2)));
        history(iter) = err;

        mu(activeSites, :) = mu_pol;

        doResidual = strcmp(stopMetric, "relres") || ...
                     (iter == 1) || ...
                     (mod(iter, residualEvery) == 0);

        if doResidual
            relres = local_matrix_free_relres(sys, problem, mu, Eext, dipoleParams);
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
                relres = local_matrix_free_relres(sys, problem, mu, Eext, dipoleParams);
                relresHistory(iter) = relres;

                if verbose
                    fprintf(' final relres = %.3e\n', relres);
                end
            end

            history = history(1:iter);
            relresHistory = relresHistory(1:iter);

            if verbose
                fprintf('SCF(iterative, matrix-free SOR) converged in %d iterations (%.2f s)\n', ...
                    iter, toc(tSCF));
            end
            break;
        end
    end

    if ~converged
        history = history(1:maxIter);
        relresHistory = relresHistory(1:maxIter);

        if isnan(relresHistory(end))
            relresHistory(end) = local_matrix_free_relres(sys, problem, mu, Eext, dipoleParams);
        end

        if verbose
            fprintf(['SCF(iterative, matrix-free SOR) hit maxIter=%d after %.2f s ' ...
                     '| final max|dmu|=%.3e\n'], ...
                maxIter, toc(tSCF), history(end));
        end
    end

    scf = struct();
    scf.method = 'iterative_sor';
    scf.variant = 'matrix_free';
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
    scf.used_thole = true;
    scf.used_matrix_solver = false;
end

function relres = local_matrix_free_relres(sys, problem, mu, Eext, dipoleParams)
    polMask = problem.polMask;
    alpha   = problem.alpha;

    Edip = thole.induced_field_from_dipoles_thole(sys, mu, dipoleParams);

    rhs = zeros(problem.nSites, 3);
    rhs(polMask, :) = alpha(polMask) .* (Eext(polMask, :) + Edip(polMask, :));

    r = rhs(polMask, :) - mu(polMask, :);

    bn = norm(rhs(polMask, :), 'fro');
    if bn == 0
        bn = 1.0;
    end

    relres = norm(r, 'fro') / bn;
end