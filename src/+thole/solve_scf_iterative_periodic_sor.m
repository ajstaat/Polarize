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
%                 .realspace_cache   optional prebuilt real-space cache
%                 .realspace_row_cache optional prebuilt row cache
%                 .kspace_cache      optional prebuilt k-space cache
%                 .kspace_mode       'auto' | 'full' | 'chunked'
%                 .kspace_memory_limit_gb positive scalar, default 8
%                 .k_block_size      positive integer, default 2048
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

    boundary = 'tinfoil';
    if isfield(ewaldParams, 'boundary') && ~isempty(ewaldParams.boundary)
        boundary = lower(char(string(ewaldParams.boundary)));
    end

    % ---------------------------------------------------------------------
    % Build or reuse caches
    
    use_thole = true;
    if isfield(scfParams, 'use_thole') && ~isempty(scfParams.use_thole)
        use_thole = logical(scfParams.use_thole);
    end
    
    rowCache = [];
    realCache = [];
    
    if isfield(scfParams, 'realspace_dipole_row_cache') && ~isempty(scfParams.realspace_dipole_row_cache)
        rowCache = scfParams.realspace_dipole_row_cache;
    elseif isfield(scfParams, 'realspace_row_cache') && ~isempty(scfParams.realspace_row_cache)
        rowCache = scfParams.realspace_row_cache;
    else
        rowCache = geom.build_periodic_realspace_dipole_row_cache(sys, problem, ewaldParams, scfParams);
    end
    
    % Legacy realCache is only needed if someone explicitly provided it
    % and you want to keep it around for diagnostics.
    if isfield(scfParams, 'realspace_cache') && ~isempty(scfParams.realspace_cache)
        realCache = scfParams.realspace_cache;
    end
    
    if isfield(scfParams, 'kspace_cache') && ~isempty(scfParams.kspace_cache)
        kCache = scfParams.kspace_cache;
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

    % ---------------------------------------------------------------------
    % Real-space row data

    row_ptr = rowCache.row_ptr;
    col_idx = rowCache.col_idx;
    dr_all = rowCache.dr;
    coeff_iso_all = rowCache.coeff_iso(:);
    coeff_dyad_all = rowCache.coeff_dyad(:);

    % ---------------------------------------------------------------------
    % Reciprocal-space metadata

    nk = kCache.num_kvec;
    kvecs = kCache.kvecs;
    pref = kCache.pref(:);
    two_pref = 2 * pref;

    % ---------------------------------------------------------------------
    % Build local diagonal blocks:
    %   Ddiag(:,:,i) = Dreal_diag(:,:,i) + Drecip_diag + Dself + Dsurf

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
        for m = 1:nk
            k = kvecs(m, :).';
            Drecip_diag = Drecip_diag + two_pref(m) * (k * k.');
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

    if nk > 0
        switch kCache.storage_mode
            case 'full'
                cos_phase = kCache.cos_phase;
                sin_phase = kCache.sin_phase;

                v = mu_pol * kvecs.';
                A = sum(cos_phase .* v, 1);
                B = sum(sin_phase .* v, 1);

            case 'chunked'
                cos_phase = [];
                sin_phase = [];
                A = zeros(1, nk);
                B = zeros(1, nk);

                pos_pol = kCache.active_pos;
                blk = kCache.k_block_size;

                for k0 = 1:blk:nk
                    k1 = min(k0 + blk - 1, nk);
                    idx = k0:k1;

                    kblk = kvecs(idx, :);
                    phase_blk = pos_pol * kblk.';
                    cos_blk = cos(phase_blk);
                    sin_blk = sin(phase_blk);
                    v_blk = mu_pol * kblk.';

                    A(idx) = sum(cos_blk .* v_blk, 1);
                    B(idx) = sum(sin_blk .* v_blk, 1);
                end

            otherwise
                error('thole:solve_scf_iterative_periodic_sor:BadKspaceMode', ...
                    'Unknown kCache.storage_mode "%s".', kCache.storage_mode);
        end
    else
        cos_phase = zeros(nPolSites, 0);
        sin_phase = zeros(nPolSites, 0);
        A = zeros(1, 0);
        B = zeros(1, 0);
    end

    history = zeros(maxIter, 1);
    relresHistory = nan(maxIter, 1);
    converged = false;

    if verbose
        if ~isempty(realCache)
            nRealPrint = realCache.nInteractions;
        else
            nRealPrint = rowCache.nInteractions;
        end
        
        fprintf(['SCF(iterative, periodic matrix-free SOR): tol=%.3e, maxIter=%d, omega=%.3f' ...
                 ' | nReal=%d | nK=%d\n'], ...
            tol, maxIter, omega, nRealPrint, kCache.num_kvec);
    end

    tSCF = tic;

    for iter = 1:maxIter
        mu_old_pol = mu_pol;

        for i = 1:nPolSites
            % -------------------------------------------------------------
            % Real-space row action for current row i

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

            % -------------------------------------------------------------
            % Reciprocal-space action for current row i

            E_recip = [0.0, 0.0, 0.0];
            if nk > 0
                switch kCache.storage_mode
                    case 'full'
                        phase_factor = cos_phase(i, :) .* A + sin_phase(i, :) .* B;
                        w = phase_factor .* two_pref.';
                        E_recip(1) = w * kvecs(:, 1);
                        E_recip(2) = w * kvecs(:, 2);
                        E_recip(3) = w * kvecs(:, 3);

                    case 'chunked'
                        ri = kCache.active_pos(i, :);
                        blk = kCache.k_block_size;

                        for k0 = 1:blk:nk
                            k1 = min(k0 + blk - 1, nk);
                            idx = k0:k1;

                            kblk = kvecs(idx, :);
                            phase_i = ri * kblk.';
                            cos_i = cos(phase_i);
                            sin_i = sin(phase_i);

                            phase_factor = cos_i .* A(idx) + sin_i .* B(idx);
                            w = phase_factor .* two_pref(idx).';

                            E_recip(1) = E_recip(1) + w * kblk(:, 1);
                            E_recip(2) = E_recip(2) + w * kblk(:, 2);
                            E_recip(3) = E_recip(3) + w * kblk(:, 3);
                        end

                    otherwise
                        error('thole:solve_scf_iterative_periodic_sor:BadKspaceMode', ...
                            'Unknown kCache.storage_mode "%s".', kCache.storage_mode);
                end
            end

            % -------------------------------------------------------------
            % Self / surface contributions

            mui_current = mu_pol(i, :);
            E_self = self_coeff * mui_current;

            E_surf = [0.0, 0.0, 0.0];
            if surf_coeff ~= 0
                Msrc = sum(mu_pol, 1);
                E_surf = surf_coeff * Msrc;
            end

            % -------------------------------------------------------------
            % Full field at row i

            E_full = Eext_pol(i, :) + E_real + E_recip + E_self + E_surf;

            Dii = Ddiag(:, :, i);
            E_off = E_full - (Dii * mui_current.').';

            rhs_i = alpha_pol(i) .* E_off;
            mu_gs_i = (Msolve(:, :, i) \ rhs_i.').';

            mu_new_i = (1 - omega) .* mui_current + omega .* mu_gs_i;

            delta_i = mu_new_i - mui_current;
            mu_pol(i, :) = mu_new_i;

            % -------------------------------------------------------------
            % Incrementally update reciprocal source sums A,B

            if nk > 0 && any(delta_i ~= 0)
                switch kCache.storage_mode
                    case 'full'
                        delta_v = delta_i * kvecs.';
                        A = A + cos_phase(i, :) .* delta_v;
                        B = B + sin_phase(i, :) .* delta_v;

                    case 'chunked'
                        ri = kCache.active_pos(i, :);
                        blk = kCache.k_block_size;

                        for k0 = 1:blk:nk
                            k1 = min(k0 + blk - 1, nk);
                            idx = k0:k1;

                            kblk = kvecs(idx, :);
                            phase_i = ri * kblk.';
                            cos_i = cos(phase_i);
                            sin_i = sin(phase_i);
                            delta_v = delta_i * kblk.';

                            A(idx) = A(idx) + cos_i .* delta_v;
                            B(idx) = B(idx) + sin_i .* delta_v;
                        end

                    otherwise
                        error('thole:solve_scf_iterative_periodic_sor:BadKspaceMode', ...
                            'Unknown kCache.storage_mode "%s".', kCache.storage_mode);
                end
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