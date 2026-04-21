function [mu, scf] = solve_scf_iterative_periodic_gmres(sys, Eext, ewaldParams, scfParams)
%SOLVE_SCF_ITERATIVE_PERIODIC_GMRES Matrix-free periodic GMRES SCF solver.
%
% [mu, scf] = thole.solve_scf_iterative_periodic_gmres(sys, Eext, ewaldParams)
% [mu, scf] = thole.solve_scf_iterative_periodic_gmres(sys, Eext, ewaldParams, scfParams)
%
% Solves the active-space linear system
%   (I - A*T) mu = A*Eext
% using MATLAB gmres with a matrix-free periodic dipole operator apply.
%
% No dense periodic interaction matrix is assembled.
%
% Inputs
%   sys          canonical polarization-system struct in atomic units
%   Eext         N x 3 external field
%   ewaldParams  struct with fields:
%                  .alpha
%                  .rcut
%                  .kcut
%                  .boundary   optional, default 'tinfoil'
%
%   scfParams    optional struct with fields:
%                  .tol
%                  .maxIter
%                  .restart
%                  .verbose
%                  .use_thole
%                  .use_preconditioner
%                  .initial_mu        optional full-space N x 3 initial guess
%                  .x0                optional N x 3 or packed 3*Npol vector
%                  .realspace_cache
%                  .realspace_row_cache
%                  .kspace_cache
%                  .kspace_mode       'auto' | 'full' | 'chunked'
%                  .kspace_memory_limit_gb
%                  .k_block_size
%
% Output
%   mu           N x 3 induced dipoles
%   scf          convergence metadata

    if nargin < 4 || isempty(scfParams)
        scfParams = struct();
    end

    io.assert_atomic_units(sys);

    if nargin < 3 || isempty(ewaldParams)
        error('thole:solve_scf_iterative_periodic_gmres:MissingEwaldParams', ...
            'ewaldParams is required.');
    end
    if ~isfield(ewaldParams, 'alpha') || isempty(ewaldParams.alpha)
        error('thole:solve_scf_iterative_periodic_gmres:MissingAlpha', ...
            'ewaldParams.alpha is required.');
    end
    if ~isfield(ewaldParams, 'rcut') || isempty(ewaldParams.rcut)
        error('thole:solve_scf_iterative_periodic_gmres:MissingRcut', ...
            'ewaldParams.rcut is required.');
    end
    if ~isfield(ewaldParams, 'kcut') || isempty(ewaldParams.kcut)
        error('thole:solve_scf_iterative_periodic_gmres:MissingKcut', ...
            'ewaldParams.kcut is required.');
    end

    problem = thole.prepare_scf_problem(sys, Eext, scfParams);

    activeSites = problem.activeSites(:);
    nPolSites   = problem.nPolSites;
    polMask     = problem.polMask;
    alpha_pol   = problem.alpha_pol(:);
    Eext_pol    = problem.Eext_pol;

    tol = 1e-10;
    if isfield(scfParams, 'tol') && ~isempty(scfParams.tol)
        tol = scfParams.tol;
    end

    maxIter = 100;
    if isfield(scfParams, 'maxIter') && ~isempty(scfParams.maxIter)
        maxIter = scfParams.maxIter;
    end

    restart = [];
    if isfield(scfParams, 'restart') && ~isempty(scfParams.restart)
        restart = scfParams.restart;
    end

    verbose = true;
    if isfield(scfParams, 'verbose') && ~isempty(scfParams.verbose)
        verbose = logical(scfParams.verbose);
    end

    use_thole = true;
    if isfield(scfParams, 'use_thole') && ~isempty(scfParams.use_thole)
        use_thole = logical(scfParams.use_thole);
    end

    use_preconditioner = true;
    if isfield(scfParams, 'use_preconditioner') && ~isempty(scfParams.use_preconditioner)
        use_preconditioner = logical(scfParams.use_preconditioner);
    end

    boundary = 'tinfoil';
    if isfield(ewaldParams, 'boundary') && ~isempty(ewaldParams.boundary)
        boundary = lower(char(string(ewaldParams.boundary)));
    end

    % ---------------------------------------------------------------------
    % Build or reuse periodic caches

    cacheParams = struct();
    cacheParams.use_thole = use_thole;
    cacheParams.verbose = verbose;

    if isfield(scfParams, 'realspace_cache') && ~isempty(scfParams.realspace_cache)
        realCache = scfParams.realspace_cache;
    else
        realCache = geom.build_periodic_realspace_cache(sys, problem, ewaldParams, cacheParams);
    end

    if isfield(scfParams, 'realspace_row_cache') && ~isempty(scfParams.realspace_row_cache)
        rowCache = scfParams.realspace_row_cache;
    else
        rowCache = geom.build_periodic_realspace_row_cache(sys, problem, realCache);
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

    dipoleParams = struct();
    dipoleParams.use_thole = use_thole;
    dipoleParams.problem = problem;
    dipoleParams.target_mask = polMask;
    dipoleParams.source_mask = polMask;
    dipoleParams.realspace_cache = realCache;
    dipoleParams.kspace_cache = kCache;

    % ---------------------------------------------------------------------
    % Right-hand side b = A * Eext in active-space vector form

    b_pol = alpha_pol .* Eext_pol;
    b = local_pack_pol(b_pol);

    % ---------------------------------------------------------------------
    % Initial guess

    if isfield(scfParams, 'x0') && ~isempty(scfParams.x0)
        x0_in = scfParams.x0;
        if isequal(size(x0_in), [problem.nSites, 3])
            x0 = local_pack_pol(x0_in(activeSites, :));
        elseif isvector(x0_in) && numel(x0_in) == 3 * nPolSites
            x0 = x0_in(:);
        else
            error('thole:solve_scf_iterative_periodic_gmres:BadX0', ...
                'scfParams.x0 must be N x 3 or a (3*Npol)x1 vector.');
        end
    else
        x0 = local_pack_pol(problem.mu0(activeSites, :));
    end

    % ---------------------------------------------------------------------
    % Optional block-diagonal preconditioner from periodic local diagonal blocks

    if use_preconditioner
        MinvBlocks = local_build_block_preconditioner(sys, problem, ewaldParams, ...
            rowCache, kCache, alpha_pol, boundary);
        M1 = @(x) local_apply_block_preconditioner(MinvBlocks, x);
    else
        M1 = [];
    end

    % ---------------------------------------------------------------------
    % Matrix-free linear operator: y = (I - A*T) x

    Afun = @(x) local_apply_linear_operator(x, sys, problem, alpha_pol, ...
        activeSites, dipoleParams, ewaldParams);

    if verbose
        fprintf(['SCF(iterative, periodic matrix-free GMRES): tol=%.3e, ' ...
                 'maxIter=%d, restart=%s, preconditioner=%d | nReal=%d | nK=%d | kMode=%s\n'], ...
            tol, maxIter, local_restart_str(restart), use_preconditioner, ...
            realCache.nInteractions, kCache.num_kvec, kCache.storage_mode);
    end

    tSolve = tic;

    if isempty(M1)
        [x, flag, relres, iter, resvec] = gmres(Afun, b, restart, tol, maxIter, [], [], x0);
    else
        [x, flag, relres, iter, resvec] = gmres(Afun, b, restart, tol, maxIter, M1, [], x0);
    end

    solveTime = toc(tSolve);

    mu_pol = local_unpack_pol(x, nPolSites);
    mu = zeros(problem.nSites, 3);
    mu(activeSites, :) = mu_pol;

    % Physical residual using the same validated periodic operator
    Edip = thole.induced_field_from_dipoles_thole_periodic(sys, mu, ewaldParams, dipoleParams);
    rhs = zeros(problem.nSites, 3);
    rhs(polMask, :) = problem.alpha(polMask) .* (Eext(polMask, :) + Edip(polMask, :));
    r = rhs(polMask, :) - mu(polMask, :);
    bn = norm(rhs(polMask, :), 'fro');
    if bn == 0
        bn = 1.0;
    end
    relres_phys = norm(r, 'fro') / bn;

    scf = struct();
    scf.method = 'gmres';
    scf.variant = 'periodic_matrix_free';
    scf.converged = (flag == 0);
    scf.flag = flag;
    scf.relres = relres;
    scf.relres_physical = relres_phys;
    scf.iter = iter;
    scf.resvec = resvec;
    scf.solve_time_s = solveTime;
    scf.used_thole = use_thole;
    scf.used_preconditioner = use_preconditioner;
    scf.nRealInteractions = realCache.nInteractions;
    scf.num_kvec = kCache.num_kvec;
    scf.kspace_storage_mode = kCache.storage_mode;

    if verbose
        fprintf(['SCF(iterative, periodic matrix-free GMRES) finished in %.2f s | ' ...
                 'flag=%d | relres(gmres)=%.3e | relres(physical)=%.3e\n'], ...
            solveTime, flag, relres, relres_phys);

        if ~isempty(resvec)
            fprintf('  residual history length = %d\n', numel(resvec));
        end
    end
end

% =========================================================================
% Local helpers
% =========================================================================

function y = local_apply_linear_operator(x, sys, problem, alpha_pol, activeSites, dipoleParams, ewaldParams)
    nPolSites = problem.nPolSites;

    mu_pol = local_unpack_pol(x, nPolSites);

    mu = zeros(problem.nSites, 3);
    mu(activeSites, :) = mu_pol;

    Edip = thole.induced_field_from_dipoles_thole_periodic(sys, mu, ewaldParams, dipoleParams);
    Edip_pol = Edip(activeSites, :);

    y_pol = mu_pol - alpha_pol .* Edip_pol;
    y = local_pack_pol(y_pol);
end

function x = local_pack_pol(mu_pol)
    x = reshape(mu_pol.', [], 1);
end

function mu_pol = local_unpack_pol(x, nPolSites)
    mu_pol = reshape(x, 3, nPolSites).';
end

function MinvBlocks = local_build_block_preconditioner(sys, problem, ewaldParams, rowCache, kCache, alpha_pol, boundary)
    nPolSites = problem.nPolSites;
    I3 = eye(3);

    % Real-space self-image diagonal blocks
    Dreal_diag = zeros(3, 3, nPolSites);

    row_ptr = rowCache.row_ptr;
    col_idx = rowCache.col_idx;
    dr_all = rowCache.dr;
    coeff_iso_all = rowCache.coeff_iso(:);
    coeff_dyad_all = rowCache.coeff_dyad(:);

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

    % Reciprocal-space diagonal block
    Drecip_diag = zeros(3, 3);
    if kCache.num_kvec > 0
        two_pref = 2 * kCache.pref(:);
        kvecs = kCache.kvecs;
        for m = 1:kCache.num_kvec
            k = kvecs(m, :).';
            Drecip_diag = Drecip_diag + two_pref(m) * (k * k.');
        end
    end

    % Self term: positive sign, matching assembled Tpol convention
    alpha_ewald = ewaldParams.alpha;
    Dself = ewald.self_tensor_block_dipole(alpha_ewald);

    % Surface term
    surf_coeff = 0.0;
    switch boundary
        case 'tinfoil'
            surf_coeff = 0.0;
        case 'vacuum'
            H = local_get_direct_lattice(sys);
            V = abs(det(H));
            surf_coeff = 4 * pi / (3 * V);
        otherwise
            error('thole:solve_scf_iterative_periodic_gmres:UnknownBoundary', ...
                'boundary must be ''tinfoil'' or ''vacuum''.');
    end
    Dsurf = surf_coeff * I3;

    MinvBlocks = zeros(3, 3, nPolSites);
    for i = 1:nPolSites
        Dii = Dreal_diag(:, :, i) + Drecip_diag + Dself + Dsurf;
        Mii = I3 - alpha_pol(i) * Dii;
        MinvBlocks(:, :, i) = inv(Mii);
    end
end

function y = local_apply_block_preconditioner(MinvBlocks, x)
    nPolSites = size(MinvBlocks, 3);
    x_pol = local_unpack_pol(x, nPolSites);
    y_pol = zeros(nPolSites, 3);

    for i = 1:nPolSites
        y_pol(i, :) = (MinvBlocks(:, :, i) * x_pol(i, :).').';
    end

    y = local_pack_pol(y_pol);
end

function s = local_restart_str(restart)
    if isempty(restart)
        s = '[]';
    else
        s = sprintf('%d', restart);
    end
end

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('Missing direct lattice on system.');
    end
end