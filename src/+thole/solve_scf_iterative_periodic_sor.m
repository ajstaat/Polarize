function [mu, scf] = solve_scf_iterative_periodic_sor(sys, Eext, ewaldParams, scfParams)
%SOLVE_SCF_ITERATIVE_PERIODIC_SOR Matrix-free periodic GS/SOR solve.
%
% [mu, scf] = thole.solve_scf_iterative_periodic_sor(sys, Eext, ewaldParams)
% [mu, scf] = thole.solve_scf_iterative_periodic_sor(sys, Eext, ewaldParams, scfParams)
%
% Solves the active-space periodic SCF equation
%   mu = alpha .* (Eext + Ereal + Erecip + Eself + Esurf)
% with a local 3x3 implicit GS/SOR sweep:
%   (I - alpha_i * D_i) mu_i = alpha_i * (Eext_i + Eoff_i)
%
% Inputs
%   sys         canonical polarization-system struct in atomic units
%   Eext        N x 3 external field
%   ewaldParams struct with fields:
%                 .alpha
%                 .rcut
%                 .kcut
%                 .boundary   optional, default 'tinfoil'
%
%   scfParams   optional solver controls:
%                 .tol
%                 .maxIter
%                 .omega
%                 .verbose
%                 .printEvery
%                 .residualEvery
%                 .stopMetric        'max_dmu' | 'relres'
%                 .use_thole         logical, default true
%                 .softening         accepted but unused in periodic solver
%                 .realspace_dipole_row_cache optional prebuilt row cache
%                 .realspace_row_cache optional prebuilt row cache
%                 .kspace_cache      optional prebuilt k-space cache
%                 .kspace_mode       'auto' | 'full' | 'blocked'
%                                    (legacy alias 'chunked' accepted)
%                 .kspace_memory_limit_gb positive scalar, default 8
%                 .k_block_size      positive integer, default 2048
%                 .use_mex_kspace    logical, default false
%
% Output
%   mu          N x 3 induced dipoles
%   scf         convergence metadata

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

    nSites      = problem.nSites;
    polMask     = problem.polMask;
    activeSites = problem.activeSites(:);
    nPolSites   = problem.nPolSites;
    alpha_pol   = problem.alpha_pol(:);
    Eext_pol    = problem.Eext_pol;

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
        error('thole:solve_scf_iterative_periodic_sor:BadOmega', ...
            'scfParams.omega must satisfy 0 < omega < 2.');
    end

    use_thole = true;
    if isfield(scfParams, 'use_thole') && ~isempty(scfParams.use_thole)
        use_thole = logical(scfParams.use_thole);
    end

    use_mex_kspace = false;
    if isfield(scfParams, 'use_mex_kspace') && ~isempty(scfParams.use_mex_kspace)
        use_mex_kspace = logical(scfParams.use_mex_kspace);
    end

    boundary = 'tinfoil';
    if isfield(ewaldParams, 'boundary') && ~isempty(ewaldParams.boundary)
        boundary = lower(char(string(ewaldParams.boundary)));
    end

    % ---------------------------------------------------------------------
    % Build or reuse real-space row cache

    if isfield(scfParams, 'realspace_dipole_row_cache') && ~isempty(scfParams.realspace_dipole_row_cache)
        rowCache = scfParams.realspace_dipole_row_cache;
    elseif isfield(scfParams, 'realspace_row_cache') && ~isempty(scfParams.realspace_row_cache)
        rowCache = scfParams.realspace_row_cache;
    else
        rowOpts = struct();
        rowOpts.profile = verbose;
        rowOpts.use_mex = false;
        rowOpts.use_thole = use_thole;
        rowCache = geom.build_active_row_cache_periodic(sys, problem, ewaldParams, rowOpts);
    end

    % ---------------------------------------------------------------------
    % Build or reuse k-space cache / plan

    if isfield(scfParams, 'kspace_cache') && ~isempty(scfParams.kspace_cache)
        kCache = scfParams.kspace_cache;
        if strcmpi(kCache.storage_mode, 'chunked')
            kCache.storage_mode = 'blocked';
        end
        if isfield(kCache, 'phase_storage') && strcmpi(kCache.phase_storage, 'chunked')
            kCache.phase_storage = 'blocked';
        end
    else
        kOpts = struct();
        if isfield(scfParams, 'kspace_mode') && ~isempty(scfParams.kspace_mode)
            kOpts.kspace_mode = scfParams.kspace_mode;
        end
        if isfield(scfParams, 'kspace_memory_limit_gb') && ~isempty(scfParams.kspace_memory_limit_gb)
            kOpts.kspace_memory_limit_gb = scfParams.kspace_memory_limit_gb;
        end
        if isfield(scfParams, 'k_block_size') && ~isempty(scfParams.k_block_size)
            kOpts.k_block_size = scfParams.k_block_size;
        end
        kOpts.verbose = verbose;
        kCache = geom.build_periodic_kspace_cache(sys, problem, ewaldParams, kOpts);
    end

    if strcmpi(kCache.storage_mode, 'chunked')
        kCache.storage_mode = 'blocked';
    end
    if ~isfield(kCache, 'two_pref') || isempty(kCache.two_pref)
        kCache.two_pref = 2 * kCache.pref(:);
    end
    if ~isfield(kCache, 'kvecs_T') || isempty(kCache.kvecs_T)
        kCache.kvecs_T = kCache.kvecs.';
    end
    if ~isfield(kCache, 'num_blocks') || isempty(kCache.num_blocks)
        [kCache.block_start, kCache.block_end, kCache.blocks] = ...
            local_make_k_blocks(kCache.kvecs, kCache.pref(:), kCache.two_pref(:), ...
                                kCache.k_block_size, kCache.active_pos, false);
        kCache.num_blocks = numel(kCache.block_start);
    end

    if use_mex_kspace && strcmpi(kCache.storage_mode, 'blocked')
        local_assert_mex_available();
    end

    % ---------------------------------------------------------------------
    % Timing instrumentation

    profileTimers = struct();
    profileTimers.t_init_recip_ab   = 0.0;
    profileTimers.t_real_rows       = 0.0;
    profileTimers.t_recip_rows      = 0.0;
    profileTimers.t_recip_updates   = 0.0;
    profileTimers.t_residual_checks = 0.0;

    % ---------------------------------------------------------------------
    % Real-space row data

    row_ptr = rowCache.row_ptr(:);
    col_idx = rowCache.col_idx(:);
    dr_all = rowCache.dr;
    coeff_iso_all = rowCache.coeff_iso(:);
    coeff_dyad_all = rowCache.coeff_dyad(:);

    % ---------------------------------------------------------------------
    % Reciprocal-space metadata

    nk = kCache.num_kvec;

    % ---------------------------------------------------------------------
    % Build local diagonal blocks

    I3 = eye(3);
    Dreal_diag = zeros(3, 3, nPolSites);

    for i = 1:nPolSites
        idx0 = row_ptr(i);
        idx1 = row_ptr(i + 1) - 1;

        if idx1 < idx0
            continue;
        end

        idx = idx0:idx1;
        selfMask = (col_idx(idx) == i);

        if any(selfMask)
            idxs = idx(selfMask);
            Dii = zeros(3, 3);

            for p = idxs
                x = dr_all(p, :).';
                Dii = Dii + coeff_iso_all(p) * I3 + coeff_dyad_all(p) * (x * x.');
            end

            Dreal_diag(:, :, i) = Dii;
        end
    end

    Drecip_diag = zeros(3, 3);
    if nk > 0
        for b = 1:kCache.num_blocks
            blk = kCache.blocks(b);
            kb = blk.kvecs;
            wp = blk.two_pref(:);

            if isempty(kb)
                continue;
            end

            for m = 1:blk.nk
                k = kb(m, :).';
                Drecip_diag = Drecip_diag + wp(m) * (k * k.');
            end
        end
    end

    alpha_ewald = ewaldParams.alpha;
    self_coeff = +(4 * alpha_ewald^3 / (3 * sqrt(pi)));
    Dself = self_coeff * I3;

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
    Dsurf = surf_coeff * I3;

    Ddiag = zeros(3, 3, nPolSites);
    Msolve = zeros(3, 3, nPolSites);
    for i = 1:nPolSites
        Ddiag(:, :, i) = Dreal_diag(:, :, i) + Drecip_diag + Dself + Dsurf;
        Msolve(:, :, i) = I3 - alpha_pol(i) * Ddiag(:, :, i);
    end

    % ---------------------------------------------------------------------
    % Initial guess

    mu = problem.mu0;
    mu_pol = mu(activeSites, :);

    % ---------------------------------------------------------------------
    % Initialize reciprocal source sums

    tTmp = tic;
    if nk > 0
        [A, B] = local_init_recip_source_sums(mu_pol, kCache);
    else
        A = zeros(1, 0);
        B = zeros(1, 0);
    end
    profileTimers.t_init_recip_ab = profileTimers.t_init_recip_ab + toc(tTmp);

    history = zeros(maxIter, 1);
    relresHistory = nan(maxIter, 1);
    converged = false;

    if verbose
        fprintf(['SCF(iterative, periodic matrix-free SOR): tol=%.3e, maxIter=%d, omega=%.3f' ...
                 ' | nReal=%d | nK=%d | kMode=%s'], ...
            tol, maxIter, omega, rowCache.nInteractions, kCache.num_kvec, kCache.storage_mode);
        if strcmpi(kCache.storage_mode, 'blocked')
            fprintf(' | mexK=%d', use_mex_kspace);
        end
        fprintf('\n');
    end

    tSCF = tic;

    for iter = 1:maxIter
        mu_old_pol = mu_pol;

        for i = 1:nPolSites
            % -------------------------------------------------------------
            % Real-space row action

            tTmp = tic;
            E_real = [0.0, 0.0, 0.0];

            idx0 = row_ptr(i);
            idx1 = row_ptr(i + 1) - 1;

            if idx1 >= idx0
                idx = idx0:idx1;
                cols = col_idx(idx);

                muNbr = mu_pol(cols, :);
                dr = dr_all(idx, :);

                muDotR = sum(muNbr .* dr, 2);
                contrib = coeff_iso_all(idx) .* muNbr + ...
                          coeff_dyad_all(idx) .* (muDotR .* dr);

                E_real = sum(contrib, 1);
            end
            profileTimers.t_real_rows = profileTimers.t_real_rows + toc(tTmp);

            % -------------------------------------------------------------
            % Reciprocal-space row action

            tTmp = tic;
            E_recip = [0.0, 0.0, 0.0];
            if nk > 0
                E_recip = local_apply_recip_row(i, A, B, kCache, use_mex_kspace);
            end
            profileTimers.t_recip_rows = profileTimers.t_recip_rows + toc(tTmp);

            % -------------------------------------------------------------
            % Self / surface

            mui_current = mu_pol(i, :);
            E_self = self_coeff * mui_current;

            E_surf = [0.0, 0.0, 0.0];
            if surf_coeff ~= 0
                Msrc = sum(mu_pol, 1);
                E_surf = surf_coeff * Msrc;
            end

            % -------------------------------------------------------------
            % Update row

            E_full = Eext_pol(i, :) + E_real + E_recip + E_self + E_surf;

            Dii = Ddiag(:, :, i);
            E_off = E_full - (Dii * mui_current.').';

            rhs_i = alpha_pol(i) .* E_off;
            mu_gs_i = (Msolve(:, :, i) \ rhs_i.').';

            mu_new_i = (1 - omega) .* mui_current + omega .* mu_gs_i;

            delta_i = mu_new_i - mui_current;
            mu_pol(i, :) = mu_new_i;

            % -------------------------------------------------------------
            % Incremental reciprocal update

            if nk > 0 && any(delta_i ~= 0)
                tTmp = tic;
                [A, B] = local_update_recip_source_sums(i, delta_i, A, B, kCache, use_mex_kspace);
                profileTimers.t_recip_updates = profileTimers.t_recip_updates + toc(tTmp);
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
            tTmp = tic;

            dipoleParams = struct();
            dipoleParams.use_thole = use_thole;
            dipoleParams.problem = problem;
            dipoleParams.target_mask = polMask;
            dipoleParams.source_mask = polMask;
            dipoleParams.realspace_row_cache = rowCache;
            dipoleParams.kspace_cache = kCache;

            Edip = thole.induced_field_from_dipoles_thole_periodic(sys, mu, ewaldParams, dipoleParams);

            rhs = zeros(nSites, 3);
            rhs(polMask, :) = problem.alpha(polMask) .* (Eext(polMask, :) + Edip(polMask, :));
            r = rhs(polMask, :) - mu(polMask, :);

            bn = norm(rhs(polMask, :), 'fro');
            if bn == 0
                bn = 1.0;
            end
            relres = norm(r, 'fro') / bn;
            relresHistory(iter) = relres;

            profileTimers.t_residual_checks = profileTimers.t_residual_checks + toc(tTmp);
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
                tTmp = tic;

                dipoleParams = struct();
                dipoleParams.use_thole = use_thole;
                dipoleParams.problem = problem;
                dipoleParams.target_mask = polMask;
                dipoleParams.source_mask = polMask;
                dipoleParams.realspace_row_cache = rowCache;
                dipoleParams.kspace_cache = kCache;

                Edip = thole.induced_field_from_dipoles_thole_periodic(sys, mu, ewaldParams, dipoleParams);
                rhs = zeros(nSites, 3);
                rhs(polMask, :) = problem.alpha(polMask) .* (Eext(polMask, :) + Edip(polMask, :));
                r = rhs(polMask, :) - mu(polMask, :);
                bn = norm(rhs(polMask, :), 'fro');
                if bn == 0
                    bn = 1.0;
                end
                relres = norm(r, 'fro') / bn;
                relresHistory(iter) = relres;

                profileTimers.t_residual_checks = profileTimers.t_residual_checks + toc(tTmp);

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
            tTmp = tic;

            dipoleParams = struct();
            dipoleParams.use_thole = use_thole;
            dipoleParams.problem = problem;
            dipoleParams.target_mask = polMask;
            dipoleParams.source_mask = polMask;
            dipoleParams.realspace_row_cache = rowCache;
            dipoleParams.kspace_cache = kCache;

            Edip = thole.induced_field_from_dipoles_thole_periodic(sys, mu, ewaldParams, dipoleParams);
            rhs = zeros(nSites, 3);
            rhs(polMask, :) = problem.alpha(polMask) .* (Eext(polMask, :) + Edip(polMask, :));
            r = rhs(polMask, :) - mu(polMask, :);
            bn = norm(rhs(polMask, :), 'fro');
            if bn == 0
                bn = 1.0;
            end
            relresHistory(end) = norm(r, 'fro') / bn;

            profileTimers.t_residual_checks = profileTimers.t_residual_checks + toc(tTmp);
        end

        if verbose
            fprintf(['SCF(iterative, periodic matrix-free SOR) hit maxIter=%d after %.2f s ' ...
                     '| final max|dmu|=%.3e\n'], ...
                maxIter, toc(tSCF), history(end));
        end
    end

    scf = struct();
    scf.method = 'iterative_periodic_sor';
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
    scf.used_thole = use_thole;
    scf.used_matrix_solver = false;
    scf.kspace_storage_mode = kCache.storage_mode;
    scf.use_mex_kspace = use_mex_kspace;
    scf.profileTimers = profileTimers;

    if verbose
        fprintf('\nPeriodic SOR timing breakdown:\n');
        fprintf('  init reciprocal A/B   : %.6f s\n', profileTimers.t_init_recip_ab);
        fprintf('  real-space row work   : %.6f s\n', profileTimers.t_real_rows);
        fprintf('  reciprocal row apply  : %.6f s\n', profileTimers.t_recip_rows);
        fprintf('  reciprocal A/B update : %.6f s\n', profileTimers.t_recip_updates);
        fprintf('  residual checks       : %.6f s\n', profileTimers.t_residual_checks);
    end
end

function [A, B] = local_init_recip_source_sums(mu_pol, kCache)
    nk = kCache.num_kvec;
    A = zeros(1, nk);
    B = zeros(1, nk);

    switch lower(kCache.storage_mode)
        case 'full'
            v = mu_pol * kCache.kvecs_T;
            A = sum(kCache.cos_phase .* v, 1);
            B = sum(kCache.sin_phase .* v, 1);

        case 'blocked'
            for b = 1:kCache.num_blocks
                blk = kCache.blocks(b);
                if blk.nk == 0
                    continue;
                end

                v_blk = mu_pol * blk.kvecs_T;
                A(blk.idx) = sum(blk.cos_phase .* v_blk, 1);
                B(blk.idx) = sum(blk.sin_phase .* v_blk, 1);
            end

        otherwise
            error('thole:solve_scf_iterative_periodic_sor:BadKspaceMode', ...
                'Unknown kCache.storage_mode "%s".', kCache.storage_mode);
    end
end

function E_recip = local_apply_recip_row(i, A, B, kCache, use_mex_kspace)
    E_recip = [0.0, 0.0, 0.0];

    switch lower(kCache.storage_mode)
        case 'full'
            phase_factor = kCache.cos_phase(i, :) .* A + kCache.sin_phase(i, :) .* B;
            w = phase_factor .* kCache.two_pref.';
            E_recip(1) = w * kCache.kvecs(:, 1);
            E_recip(2) = w * kCache.kvecs(:, 2);
            E_recip(3) = w * kCache.kvecs(:, 3);

        case 'blocked'
            for b = 1:kCache.num_blocks
                blk = kCache.blocks(b);
                if blk.nk == 0
                    continue;
                end

                if use_mex_kspace
                    Eblk = mex_periodic_kspace_block( ...
                        'apply_row', ...
                        A(blk.idx), ...
                        B(blk.idx), ...
                        blk.cos_phase(i, :), ...
                        blk.sin_phase(i, :), ...
                        blk.two_pref, ...
                        blk.kvecs);
                    E_recip = E_recip + Eblk;
                else
                    cos_i = blk.cos_phase(i, :);
                    sin_i = blk.sin_phase(i, :);

                    phase_factor = cos_i .* A(blk.idx) + sin_i .* B(blk.idx);
                    w = phase_factor .* blk.two_pref.';

                    E_recip(1) = E_recip(1) + w * blk.kvecs(:, 1);
                    E_recip(2) = E_recip(2) + w * blk.kvecs(:, 2);
                    E_recip(3) = E_recip(3) + w * blk.kvecs(:, 3);
                end
            end

        otherwise
            error('thole:solve_scf_iterative_periodic_sor:BadKspaceMode', ...
                'Unknown kCache.storage_mode "%s".', kCache.storage_mode);
    end
end

function [A, B] = local_update_recip_source_sums(i, delta_i, A, B, kCache, use_mex_kspace)
    switch lower(kCache.storage_mode)
        case 'full'
            delta_v = delta_i * kCache.kvecs_T;
            A = A + kCache.cos_phase(i, :) .* delta_v;
            B = B + kCache.sin_phase(i, :) .* delta_v;

        case 'blocked'
            for b = 1:kCache.num_blocks
                blk = kCache.blocks(b);
                if blk.nk == 0
                    continue;
                end

                if use_mex_kspace
                    [Ablk, Bblk] = mex_periodic_kspace_block( ...
                        'update_ab', ...
                        A(blk.idx), ...
                        B(blk.idx), ...
                        delta_i, ...
                        blk.cos_phase(i, :), ...
                        blk.sin_phase(i, :), ...
                        blk.kvecs_T);
                    A(blk.idx) = Ablk;
                    B(blk.idx) = Bblk;
                else
                    cos_i = blk.cos_phase(i, :);
                    sin_i = blk.sin_phase(i, :);
                    delta_v = delta_i * blk.kvecs_T;

                    A(blk.idx) = A(blk.idx) + cos_i .* delta_v;
                    B(blk.idx) = B(blk.idx) + sin_i .* delta_v;
                end
            end

        otherwise
            error('thole:solve_scf_iterative_periodic_sor:BadKspaceMode', ...
                'Unknown kCache.storage_mode "%s".', kCache.storage_mode);
    end
end

function [blockStart, blockEnd, blocks] = local_make_k_blocks(kvecs, pref, two_pref, k_block_size, pos_pol, store_phase)
    nk = size(kvecs, 1);
    nPolSites = size(pos_pol, 1);

    if nk == 0
        blockStart = zeros(0, 1);
        blockEnd = zeros(0, 1);
        blocks = repmat(local_empty_block(nPolSites), 0, 1);
        return;
    end

    blockStart = (1:k_block_size:nk).';
    nBlocks = numel(blockStart);
    blockEnd = zeros(nBlocks, 1);
    blocks = repmat(local_empty_block(nPolSites), nBlocks, 1);

    for b = 1:nBlocks
        i0 = blockStart(b);
        i1 = min(i0 + k_block_size - 1, nk);
        idx = i0:i1;

        blockEnd(b) = i1;
        blocks(b).idx = idx;
        blocks(b).kvecs = kvecs(idx, :);
        blocks(b).kvecs_T = kvecs(idx, :).';
        blocks(b).pref = pref(idx);
        blocks(b).two_pref = two_pref(idx);
        blocks(b).nk = numel(idx);

        if store_phase
            phase = pos_pol * blocks(b).kvecs_T;
            blocks(b).cos_phase = cos(phase);
            blocks(b).sin_phase = sin(phase);
        else
            blocks(b).cos_phase = zeros(nPolSites, 0);
            blocks(b).sin_phase = zeros(nPolSites, 0);
        end
    end
end

function blk = local_empty_block(nPolSites)
    if nargin < 1
        nPolSites = 0;
    end
    blk = struct( ...
        'idx', zeros(1, 0), ...
        'kvecs', zeros(0, 3), ...
        'kvecs_T', zeros(3, 0), ...
        'pref', zeros(0, 1), ...
        'two_pref', zeros(0, 1), ...
        'nk', 0, ...
        'cos_phase', zeros(nPolSites, 0), ...
        'sin_phase', zeros(nPolSites, 0));
end

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('thole:solve_scf_iterative_periodic_sor:MissingLattice', ...
            'Missing direct lattice on system.');
    end
end

function local_assert_mex_available()
    if exist('mex_periodic_kspace_block', 'file') ~= 3
        error('thole:solve_scf_iterative_periodic_sor:MissingMexKspace', ...
            ['use_mex_kspace=true but mex_periodic_kspace_block is not available.\n' ...
             'Compile it first, e.g.:\n' ...
             '  mex -O src/+thole/private/mex_periodic_kspace_block.c -outdir src/+thole/private']);
    end
end