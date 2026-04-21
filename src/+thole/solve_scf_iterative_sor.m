function [mu, scf] = solve_scf_iterative_sor(sys, Eext, scfParams)
%SOLVE_SCF_ITERATIVE_SOR Matrix-free GS/SOR solve using row-wise cached neighbors.
%
% This solver does NOT assemble dense Tpol. It performs row sweeps over a
% directed active-space neighbor cache.
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
%   .row_cache         optional prebuilt active row cache
%   .geom_cache        optional legacy unordered pair cache
%   .checkResidualAgainstLegacy   logical, default false
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

    checkResidualAgainstLegacy = false;
    if isfield(scfParams, 'checkResidualAgainstLegacy') && ~isempty(scfParams.checkResidualAgainstLegacy)
        checkResidualAgainstLegacy = logical(scfParams.checkResidualAgainstLegacy);
    end

    % ---------------------------------------------------------------------
    % Row cache for actual solver work
    % Priority:
    %   1) caller-provided row_cache
    %   2) caller-provided legacy geom_cache -> convert to row cache
    %   3) build row cache directly
    % ---------------------------------------------------------------------
    useCallerRowCache = isfield(scfParams, 'row_cache') && ~isempty(scfParams.row_cache);
    useCallerGeomCache = isfield(scfParams, 'geom_cache') && ~isempty(scfParams.geom_cache);

    if useCallerRowCache
        rowCache = scfParams.row_cache;
    elseif useCallerGeomCache
        pairCache = scfParams.geom_cache;

        requestedRcut = inf;
        if isfield(scfParams, 'rcut') && ~isempty(scfParams.rcut)
            requestedRcut = scfParams.rcut;
        end

        if isfield(pairCache, 'rcut') && isfinite(requestedRcut) && pairCache.rcut + 1e-12 < requestedRcut
            error('thole:solve_scf_iterative_sor:CacheRcutMismatch', ...
                'Provided geom_cache has rcut=%.6g but solver requested rcut=%.6g.', ...
                pairCache.rcut, requestedRcut);
        end

        rowCache = geom.build_active_row_cache(sys, problem, pairCache);
    else
        rowOpts = struct();
        if isfield(scfParams, 'rcut') && ~isempty(scfParams.rcut)
            rowOpts.rcut = scfParams.rcut;
        else
            rowOpts.rcut = inf;
        end

        rowCache = geom.build_active_row_cache(sys, problem, rowOpts);
    end

    % ---------------------------------------------------------------------
    % Legacy unordered pair cache only for optional residual cross-check
    % ---------------------------------------------------------------------
    if checkResidualAgainstLegacy
        if useCallerGeomCache
            pairCacheLegacy = scfParams.geom_cache;
        else
            cacheOpts = struct();
            cacheOpts.site_mask = polMask;

            if isfield(scfParams, 'rcut') && ~isempty(scfParams.rcut)
                cacheOpts.rcut = scfParams.rcut;
            else
                cacheOpts.rcut = inf;
            end

            pairCacheLegacy = geom.build_nonperiodic_pair_cache(sys, cacheOpts);
        end

        dipoleParams = struct();
        dipoleParams.exclude_self = true;
        dipoleParams.softening = softening;
        dipoleParams.target_mask = polMask;
        dipoleParams.source_mask = polMask;
        if isfield(scfParams, 'rcut') && ~isempty(scfParams.rcut)
            dipoleParams.rcut = scfParams.rcut;
        end
        dipoleParams.geom_cache = pairCacheLegacy;
    else
        dipoleParams = struct(); %#ok<NASGU>
    end

    mu = problem.mu0;
    mu_pol = mu(activeSites, :);

    history = zeros(maxIter, 1);
    relresHistory = nan(maxIter, 1);
    converged = false;

    row_ptr = rowCache.row_ptr;
    col_idx = rowCache.col_idx;
    dr_all  = rowCache.dr;

    haveThole = isfield(rowCache, 'thole_f3') && isfield(rowCache, 'thole_f5');
    if haveThole
        f3_all = rowCache.thole_f3;
        f5_all = rowCache.thole_f5;
    else
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

                muNbr = mu_pol(cols, :);
                dr    = dr_all(idx, :);

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

                coeff1 = 3 .* (f5 .* muDotR .* invR5);
                coeff2 = (f3 .* invR3);

                contrib = coeff1 .* dr - coeff2 .* muNbr;
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
            relres = local_matrix_free_relres_rowcache(problem, mu_pol, Eext_pol, rowCache, softening, sys);

            if checkResidualAgainstLegacy
                relresLegacy = local_matrix_free_relres_legacy(sys, problem, mu, Eext, dipoleParams);
                relDiff = abs(relres - relresLegacy) / max(1.0, abs(relresLegacy));

                if relDiff > 1e-11
                    error('thole:solve_scf_iterative_sor:ResidualMismatch', ...
                        ['Row-cache residual and legacy residual disagree. ' ...
                         'row=%.16e legacy=%.16e relDiff=%.3e'], ...
                        relres, relresLegacy, relDiff);
                end
            end

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
                relres = local_matrix_free_relres_rowcache(problem, mu_pol, Eext_pol, rowCache, softening, sys);
                relresHistory(iter) = relres;

                if checkResidualAgainstLegacy
                    relresLegacy = local_matrix_free_relres_legacy(sys, problem, mu, Eext, dipoleParams);
                    relDiff = abs(relres - relresLegacy) / max(1.0, abs(relresLegacy));
                    if relDiff > 1e-11
                        error('thole:solve_scf_iterative_sor:ResidualMismatch', ...
                            ['Row-cache residual and legacy residual disagree. ' ...
                             'row=%.16e legacy=%.16e relDiff=%.3e'], ...
                            relres, relresLegacy, relDiff);
                    end
                end

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
            relresHistory(end) = local_matrix_free_relres_rowcache(problem, mu_pol, Eext_pol, rowCache, softening, sys);

            if checkResidualAgainstLegacy
                relresLegacy = local_matrix_free_relres_legacy(sys, problem, mu, Eext, dipoleParams);
                relDiff = abs(relresHistory(end) - relresLegacy) / max(1.0, abs(relresLegacy));
                if relDiff > 1e-11
                    error('thole:solve_scf_iterative_sor:ResidualMismatch', ...
                        ['Row-cache residual and legacy residual disagree. ' ...
                         'row=%.16e legacy=%.16e relDiff=%.3e'], ...
                        relresHistory(end), relresLegacy, relDiff);
                end
            end
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

function relres = local_matrix_free_relres_rowcache(problem, mu_pol, Eext_pol, rowCache, softening, sys)
    alpha_pol = problem.alpha_pol(:);
    nPolSites = problem.nPolSites;
    activeSites = problem.activeSites(:);

    row_ptr = rowCache.row_ptr;
    col_idx = rowCache.col_idx;
    dr_all  = rowCache.dr;

    haveThole = isfield(rowCache, 'thole_f3') && isfield(rowCache, 'thole_f5');
    if haveThole
        f3_all = rowCache.thole_f3;
        f5_all = rowCache.thole_f5;
    else
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

    Edip_pol = zeros(nPolSites, 3);

    for i = 1:nPolSites
        k0 = row_ptr(i);
        k1 = row_ptr(i + 1) - 1;

        if k1 < k0
            continue;
        end

        idx = k0:k1;
        cols = col_idx(idx);

        muNbr = mu_pol(cols, :);
        dr    = dr_all(idx, :);

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

        coeff1 = 3 .* (f5 .* muDotR .* invR5);
        coeff2 = (f3 .* invR3);

        contrib = coeff1 .* dr - coeff2 .* muNbr;
        Edip_pol(i, :) = sum(contrib, 1);
    end

    rhs_pol = alpha_pol .* (Eext_pol + Edip_pol);
    r_pol = rhs_pol - mu_pol;

    bn = norm(rhs_pol, 'fro');
    if bn == 0
        bn = 1.0;
    end

    relres = norm(r_pol, 'fro') / bn;
end

function relres = local_matrix_free_relres_legacy(sys, problem, mu, Eext, dipoleParams)
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