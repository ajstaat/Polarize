function [mu, scf] = solve_scf_sor(sys, Eext, T, scfParams)
%SOLVE_SCF_SOR Solve (I - A*T) mu = A*Eext by block GS / SOR sweeps.
%
% This implementation is derived directly from the block rows of
%   M * mu = b,  where  M = I - A*T,  b = A*Eext
%
% Vector ordering matches util.stack_xyz / util.unstack_xyz:
%   [x1; y1; z1; x2; y2; z2; ...]
%
% omega = 1.0   -> Gauss-Seidel
% 0 < omega < 1 -> under-relaxed GS
% 1 < omega < 2 -> SOR

    nSites = sys.n_sites;

    verbose = false;
    if isfield(scfParams, 'verbose') && ~isempty(scfParams.verbose)
        verbose = logical(scfParams.verbose);
    end

    printEvery = 10;
    if isfield(scfParams, 'printEvery') && ~isempty(scfParams.printEvery)
        printEvery = scfParams.printEvery;
    end

    if ~isequal(size(Eext), [nSites, 3])
        error('Eext must be N x 3.');
    end
    if ~isequal(size(T), [3*nSites, 3*nSites])
        error('T must be 3N x 3N.');
    end

    tol = 1e-8;
    if isfield(scfParams, 'tol') && ~isempty(scfParams.tol)
        tol = scfParams.tol;
    end

    maxIter = 1000;
    if isfield(scfParams, 'maxIter') && ~isempty(scfParams.maxIter)
        maxIter = scfParams.maxIter;
    end

    omega = 1.0;
    if isfield(scfParams, 'omega') && ~isempty(scfParams.omega)
        omega = scfParams.omega;
    end

    if omega <= 0 || omega >= 2
        error('scfParams.omega must satisfy 0 < omega < 2.');
    end

    % Build the exact same linear system used by the direct solver:
    %   M = I - A*T
    %   b = A*Eext_vec
    A = thole.build_alpha_matrix(sys);
    Eext_vec = util.stack_xyz(Eext);
    M = eye(3*nSites) - A*T;
    b = A * Eext_vec;

    if isfield(scfParams, 'initial_mu') && ~isempty(scfParams.initial_mu)
        mu0 = scfParams.initial_mu;
        if ~isequal(size(mu0), [nSites, 3])
            error('scfParams.initial_mu must be N x 3.');
        end
        mu_vec = util.stack_xyz(mu0);
    else
        mu_vec = zeros(3*nSites, 1);
    end

    bn = norm(b);
    if bn == 0
        bn = 1.0;
    end

    history = zeros(maxIter, 1);
    converged = false;

    % Precompute 3x3 block indices.
    activeIdx = find(logical(sys.site_is_polarizable(:)));
    blocks = cell(numel(activeIdx), 1);
    for k = 1:numel(activeIdx)
        i = activeIdx(k);
        blocks{k} = (3*i-2):(3*i);
    end

    for iter = 1:maxIter
        mu_old = mu_vec;

        % Block GS / SOR sweep:
        % M_ii * mu_i^(GS) =
        %   b_i
        %   - sum_{j<i} M_ij * mu_j^(new)
        %   - sum_{j>i} M_ij * mu_j^(old)
        for k = 1:numel(activeIdx)
            Ii = blocks{k};

            rhs_i = b(Ii);

            % Earlier blocks: use newest values from mu_vec
            for j = 1:(k-1)
                Ij = blocks{j};
                rhs_i = rhs_i - M(Ii, Ij) * mu_vec(Ij);
            end

            % Later blocks: use old values from mu_old
            for j = (k+1):numel(activeIdx)
                Ij = blocks{j};
                rhs_i = rhs_i - M(Ii, Ij) * mu_old(Ij);
            end

            mu_gs_i = M(Ii, Ii) \ rhs_i;
            mu_vec(Ii) = (1 - omega) * mu_old(Ii) + omega * mu_gs_i;
        end

        r = b - M * mu_vec;
        relres = norm(r) / bn;
        history(iter) = relres;

        if verbose && (iter == 1 || mod(iter, printEvery) == 0 || relres < tol)
            fprintf('  iter %4d | relres = %.3e\n', iter, relres);
        end

        if relres < tol
            converged = true;
            history = history(1:iter);
            break;
        end
    end

    if verbose
        if converged
            fprintf('SCF(SOR) converged in %d iterations | relres = %.3e\n', ...
                numel(history), history(end));
        else
            fprintf('SCF(SOR) hit maxIter=%d | final relres = %.3e\n', ...
                maxIter, history(end));
        end
    end

    if ~converged
        history = history(1:maxIter);
    end

    mu = util.unstack_xyz(mu_vec);

    scf = struct();
    scf.method = 'sor';
    scf.omega = omega;
    scf.converged = converged;
    scf.nIter = numel(history);
    scf.relres = history(end);
    scf.history = history;
    scf.used_matrix_solver = true;
end