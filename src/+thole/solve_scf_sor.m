function [mu, scf] = solve_scf_sor(sys, Eext, T, scfParams)
%SOLVE_SCF_SOR Solve induced dipoles by block GS / SOR sweeps on T.
%
% Uses the fixed-point form
%   mu_i = alpha_i * (Eext_i + sum_j T_ij mu_j)
%
% with Gauss-Seidel / SOR ordering:
%   earlier j use newest values
%   later   j use previous-iteration values
%
% omega = 1.0 -> GS
% 0 < omega < 1 -> under-relaxed GS
% 1 < omega < 2 -> SOR

    nSites = sys.n_sites;

    if ~isequal(size(Eext), [nSites, 3])
        error('Eext must be N x 3.');
    end
    if ~isequal(size(T), [3*nSites, 3*nSites])
        error('T must be 3N x 3N.');
    end

    % ---------------------------
    % Parse options
    % ---------------------------
    tol = 1e-8;
    if isfield(scfParams,'tol') && ~isempty(scfParams.tol)
        tol = scfParams.tol;
    end

    maxIter = 1000;
    if isfield(scfParams,'maxIter') && ~isempty(scfParams.maxIter)
        maxIter = scfParams.maxIter;
    end

    omega = 1.0;
    if isfield(scfParams,'omega') && ~isempty(scfParams.omega)
        omega = scfParams.omega;
    end
    if omega <= 0 || omega >= 2
        error('scfParams.omega must satisfy 0 < omega < 2.');
    end

    verbose = false;
    if isfield(scfParams,'verbose') && ~isempty(scfParams.verbose)
        verbose = logical(scfParams.verbose);
    end

    printEvery = 10;
    if isfield(scfParams,'printEvery') && ~isempty(scfParams.printEvery)
        printEvery = scfParams.printEvery;
    end

    residualEvery = 25;
    if isfield(scfParams,'residualEvery') && ~isempty(scfParams.residualEvery)
        residualEvery = scfParams.residualEvery;
    end

    polMask = logical(sys.site_is_polarizable(:));
    alpha   = sys.site_alpha(:);

    if numel(polMask) ~= nSites || numel(alpha) ~= nSites
        error('sys.site_is_polarizable and sys.site_alpha must both have length nSites.');
    end

    if isfield(scfParams,'initial_mu') && ~isempty(scfParams.initial_mu)
        mu = scfParams.initial_mu;
        if ~isequal(size(mu), [nSites, 3])
            error('scfParams.initial_mu must be N x 3.');
        end
    else
        mu = zeros(nSites, 3);
    end

    activeIdx = find(polMask);
    nActive = numel(activeIdx);

    % Precompute block indices once.
    blkStart = 3*activeIdx - 2;
    blkMid   = 3*activeIdx - 1;
    blkEnd   = 3*activeIdx;
    blocks   = [blkStart, blkMid, blkEnd];

    % Precompute right-hand normalization vector b = A*Eext for residuals.
    % Only used for occasional diagnostics.
    b = zeros(3*nSites, 1);
    for k = 1:nActive
        i = activeIdx(k);
        Ii = blocks(k,1):blocks(k,3);
        b(Ii) = alpha(i) * Eext(i,:).';
    end
    bn = norm(b);
    if bn == 0
        bn = 1.0;
    end

    history = zeros(maxIter,1);
    relresHistory = nan(maxIter,1);
    converged = false;

    if verbose
        if abs(omega - 1.0) < 1e-14
            solverName = 'GS';
        else
            solverName = 'SOR';
        end
        fprintf('SCF(%s): tol=%.3e, maxIter=%d, omega=%.3f\n', ...
            solverName, tol, maxIter, omega);
    end

    for iter = 1:maxIter
        mu_old = mu;

        for k = 1:nActive
            i = activeIdx(k);
            Ii = blocks(k,1):blocks(k,3);

            % Start from external field at site i
            E_loc = Eext(i,:).';

            % Earlier active sites: use newest mu
            for j = 1:(k-1)
                jj = activeIdx(j);
                Jj = blocks(j,1):blocks(j,3);
                E_loc = E_loc + T(Ii, Jj) * mu(jj,:).';
            end

            % Later active sites: use previous-iteration mu
            for j = (k+1):nActive
                jj = activeIdx(j);
                Jj = blocks(j,1):blocks(j,3);
                E_loc = E_loc + T(Ii, Jj) * mu_old(jj,:).';
            end

            mu_gs_i = alpha(i) * E_loc;
            mu(i,:) = ((1 - omega) * mu_old(i,:).' + omega * mu_gs_i).';
        end

        dmu = mu - mu_old;
        err = max(sqrt(sum(dmu.^2, 2)));
        history(iter) = err;

        doResidual = (iter == 1) || (mod(iter, residualEvery) == 0) || (err < tol);
        if doResidual
            relres = local_relres_from_T(mu, T, alpha, polMask, Eext, b, bn);
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

        if err < tol
            converged = true;
            history = history(1:iter);
            relresHistory = relresHistory(1:iter);
            break;
        end
    end

    if ~converged
        history = history(1:maxIter);
        relresHistory = relresHistory(1:maxIter);
    end

    scf = struct();
    scf.method = 'sor';
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
    scf.used_matrix_solver = true;
end

function relres = local_relres_from_T(mu, T, alpha, polMask, Eext, b, bn)
% Compute || b - (I - A*T) mu || / ||b|| without building A or M.

    nSites = size(mu,1);
    mu_vec = util.stack_xyz(mu);
    Tmu = T * mu_vec;

    r = b - mu_vec;

    activeIdx = find(polMask);
    for k = 1:numel(activeIdx)
        i = activeIdx(k);
        Ii = (3*i-2):(3*i);
        r(Ii) = r(Ii) + alpha(i) * Tmu(Ii);
    end

    relres = norm(r) / bn;
end