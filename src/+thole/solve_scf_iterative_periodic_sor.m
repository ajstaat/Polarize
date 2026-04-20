function [mu, scf] = solve_scf_iterative_periodic_sor(sys, Eext, ewaldParams, scfParams)
%SOLVE_SCF_ITERATIVE_PERIODIC_SOR Matrix-free periodic GS/SOR SCF solver.
%
% [mu, scf] = thole.solve_scf_iterative_periodic_sor(sys, Eext, ewaldParams)
% [mu, scf] = thole.solve_scf_iterative_periodic_sor(sys, Eext, ewaldParams, scfParams)
%
% This solver does NOT assemble the periodic dipole interaction matrix.
% It performs row sweeps over:
%   - periodic real-space row cache
%   - incremental reciprocal-space structure-factor sums
%   - analytic self term
%   - optional surface term
%
% Crucially, this version treats the local periodic diagonal block
% implicitly via a 3x3 solve per site:
%
%   (I - alpha_i * D_i) mu_i = alpha_i (Eext_i + Eoff_i)
%
% where D_i is the on-site operator block coming from:
%   - real-space self-image terms
%   - reciprocal-space diagonal block
%   - analytic self block
%   - surface self block
%
% This is the key fix over the earlier GS/SOR-like explicit version.

    if nargin < 4 || isempty(scfParams)
        scfParams = struct();
    end

    io.assert_atomic_units(sys);

    if nargin < 3 || isempty(ewaldParams)
        error('thole:solve_scf_iterative_periodic_sor:MissingEwaldParams', ...
            'ewaldParams is required.');
    end
    if ~isfield(ewaldParams, 'alpha') || isempty(ewaldParams.alpha)
        error('thole:solve_scf_iterative_periodic_sor:MissingAlpha', ...
            'ewaldParams.alpha is required.');
    end
    if ~isfield(ewaldParams, 'rcut') || isempty(ewaldParams.rcut)
        error('thole:solve_scf_iterative_periodic_sor:MissingRcut', ...
            'ewaldParams.rcut is required.');
    end
    if ~isfield(ewaldParams, 'kcut') || isempty(ewaldParams.kcut)
        error('thole:solve_scf_iterative_periodic_sor:MissingKcut', ...
            'ewaldParams.kcut is required.');
    end

    problem = thole.prepare_scf_problem(sys, Eext, scfParams);

    nSites        = problem.nSites;
    polMask       = problem.polMask;
    activeSites   = problem.activeSites(:);
    nPolSites     = problem.nPolSites;
    alpha_pol     = problem.alpha_pol(:);
    Eext_pol      = problem.Eext_pol;

    tol           = problem.tol;
    maxIter       = problem.maxIter;
    verbose       = problem.verbose;
    printEvery    = problem.printEvery;
    residualEvery = problem.residualEvery;

    omega = 1.0;
    if isfield(scfParams, 'omega') && ~isempty(scfParams.omega)
        omega = scfParams.omega;
    end
    if omega <= 0 || omega >= 2
        error('thole:solve_scf_iterative_periodic_sor:BadOmega', ...
            'scfParams.omega must satisfy 0 < omega < 2.');
    end

    stopMetric = "max_dmu";
    if isfield(problem, 'stopMetric') && ~isempty(problem.stopMetric)
        stopMetric = lower(string(problem.stopMetric));
    end

    use_thole = true;
    if isfield(scfParams, 'use_thole') && ~isempty(scfParams.use_thole)
        use_thole = logical(scfParams.use_thole);
    end

    boundary = 'tinfoil';
    if isfield(ewaldParams, 'boundary') && ~isempty(ewaldParams.boundary)
        boundary = lower(char(string(ewaldParams.boundary)));
    end

    % ---------------------------------------------------------------------
    % Build caches once

    cacheParams = struct();
    cacheParams.use_thole = use_thole;

    realCache = geom.build_periodic_realspace_cache(sys, problem, ewaldParams, cacheParams);
    rowCache  = geom.build_periodic_realspace_row_cache(sys, problem, realCache);
    kCache    = geom.build_periodic_kspace_cache(sys, problem, ewaldParams);

    mu = problem.mu0;
    mu_pol = mu(activeSites, :);

    % ---------------------------------------------------------------------
    % Reciprocal-space incremental sums

    nk = kCache.num_kvec;
    if nk > 0
        kvecs = kCache.kvecs;                % Nk x 3
        cos_phase = kCache.cos_phase;        % nPol x Nk
        sin_phase = kCache.sin_phase;        % nPol x Nk
        pref = kCache.pref(:);               % Nk x 1
        two_pref = 2 * pref;                 % Nk x 1

        v = mu_pol * kvecs.';                % nPol x Nk
        A = sum(cos_phase .* v, 1);          % 1 x Nk
        B = sum(sin_phase .* v, 1);          % 1 x Nk
    else
        kvecs = zeros(0,3);
        cos_phase = zeros(nPolSites,0);
        sin_phase = zeros(nPolSites,0);
        two_pref = zeros(0,1);
        A = zeros(1,0);
        B = zeros(1,0);
    end

    % Surface-term source sum over active polarizable sites
    Msrc = sum(mu_pol, 1);

    alpha_ewald = ewaldParams.alpha;
    self_coeff = -(4 * alpha_ewald^3 / (3 * sqrt(pi)));

    surf_coeff = 0.0;
    switch boundary
        case 'tinfoil'
            surf_coeff = 0.0;
        case 'vacuum'
            H = local_get_direct_lattice(sys);
            V = abs(det(H));
            surf_coeff = 4 * pi / (3 * V);
        otherwise
            error('thole:solve_scf_iterative_periodic_sor:UnknownBoundary', ...
                'boundary must be ''tinfoil'' or ''vacuum''.');
    end

    % ---------------------------------------------------------------------
    % Precompute local diagonal blocks D_i and local solves
    %
    % D_i = D_real_selfimage_i + D_recip_diag + D_self + D_surf
    %
    % Then:
    %   (I - alpha_i * D_i) mu_i = alpha_i * (Eext_i + Eoff_i)

    I3 = eye(3);

    % Real-space self-image diagonal blocks: site-dependent because Thole
    % damping can depend on alpha_i.
    Dreal_diag = zeros(3, 3, nPolSites);

    row_ptr = rowCache.row_ptr;
    col_idx = rowCache.col_idx;
    dr_all = rowCache.dr;
    coeff_iso_all = rowCache.coeff_iso(:);
    coeff_dyad_all = rowCache.coeff_dyad(:);

    for i = 1:nPolSites
        k0 = row_ptr(i);
        k1 = row_ptr(i + 1) - 1;

        if k1 < k0
            continue;
        end

        idx = k0:k1;
        selfMask = (col_idx(idx) == i);

        if any(selfMask)
            idxs = idx(selfMask);
            Dii = zeros(3,3);

            for p = idxs
                x = dr_all(p, :).';
                Dii = Dii + coeff_iso_all(p) * I3 + coeff_dyad_all(p) * (x * x.');
            end

            Dreal_diag(:,:,i) = Dii;
        end
    end

    % Reciprocal-space diagonal block is site-independent because
    % cos(k·(r_i-r_i)) = 1 for every site.
    Drecip_diag = zeros(3,3);
    if nk > 0
        for m = 1:nk
            k = kvecs(m, :).';
            Drecip_diag = Drecip_diag + two_pref(m) * (k * k.');
        end
    end

    Dself = self_coeff * I3;
    Dsurf = surf_coeff * I3;

    Ddiag = zeros(3,3,nPolSites);
    Msolve = zeros(3,3,nPolSites);

    for i = 1:nPolSites
        Ddiag(:,:,i) = Dreal_diag(:,:,i) + Drecip_diag + Dself + Dsurf;
        Msolve(:,:,i) = I3 - alpha_pol(i) * Ddiag(:,:,i);
    end

    % ---------------------------------------------------------------------
    % Histories / residual setup

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
        fprintf(['SCF(iterative, periodic matrix-free SOR): tol=%.3e, maxIter=%d, ' ...
                 'omega=%.3f | nReal=%d | nK=%d\n'], ...
            tol, maxIter, omega, realCache.nInteractions, kCache.num_kvec);
    end

    tSCF = tic;

    for iter = 1:maxIter
        mu_old_pol = mu_pol;

        for i = 1:nPolSites
            % -------------------------------------------------------------
            % Full current row field components

            E_real = [0.0, 0.0, 0.0];

            k0 = row_ptr(i);
            k1 = row_ptr(i + 1) - 1;

            if k1 >= k0
                idx = k0:k1;
                cols = col_idx(idx);
                muNbr = mu_pol(cols, :);
                dr = dr_all(idx, :);

                muDotR = sum(muNbr .* dr, 2);
                contrib = coeff_iso_all(idx) .* muNbr + ...
                          coeff_dyad_all(idx) .* (muDotR .* dr);

                E_real = sum(contrib, 1);
            end

            E_recip = [0.0, 0.0, 0.0];
            if nk > 0
                phase_factor = cos_phase(i, :) .* A + sin_phase(i, :) .* B;   % 1 x Nk
                w = phase_factor .* two_pref.';                               % 1 x Nk

                E_recip(1) = w * kvecs(:,1);
                E_recip(2) = w * kvecs(:,2);
                E_recip(3) = w * kvecs(:,3);
            end

            mui_current = mu_pol(i, :);

            E_self = self_coeff * mui_current;
            E_surf = surf_coeff * Msrc;

            E_full = Eext_pol(i, :) + E_real + E_recip + E_self + E_surf;

            % -------------------------------------------------------------
            % Remove local diagonal action and solve local 3x3 block

            Dii = Ddiag(:,:,i);
            E_off = E_full - (Dii * mui_current.').';

            rhs = alpha_pol(i) * E_off;
            mu_gs = (Msolve(:,:,i) \ rhs.').';

            % Relaxed SOR update
            mu_new_i = (1 - omega) .* mui_current + omega .* mu_gs;
            delta_i = mu_new_i - mui_current;

            if any(delta_i ~= 0)
                mu_pol(i, :) = mu_new_i;

                % Update reciprocal-space source sums incrementally
                if nk > 0
                    delta_v = delta_i * kvecs.';          % 1 x Nk
                    A = A + cos_phase(i, :) .* delta_v;
                    B = B + sin_phase(i, :) .* delta_v;
                end

                % Update surface source sum incrementally
                Msrc = Msrc + delta_i;
            end
        end

        dmu_pol = mu_pol - mu_old_pol;
        err = max(sqrt(sum(dmu_pol.^2, 2)));
        history(iter) = err;

        mu(activeSites, :) = mu_pol;

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
                fprintf('SCF(iterative, periodic matrix-free SOR) converged in %d iterations (%.2f s)\n', ...
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
            fprintf(['SCF(iterative, periodic matrix-free SOR) hit maxIter=%d after %.2f s ' ...
                     '| final max|dmu| = %.3e\n'], ...
                maxIter, toc(tSCF), history(end));
        end
    end

    scf = struct();
    scf.method = 'iterative_sor';
    scf.variant = 'periodic_matrix_free';
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
    scf.used_thole = use_thole;
    scf.used_matrix_solver = false;
    scf.nRealInteractions = realCache.nInteractions;
    scf.num_kvec = kCache.num_kvec;
end

function relres = local_periodic_matrix_free_relres(sys, problem, mu, Eext, ewaldParams, dipoleParams)
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