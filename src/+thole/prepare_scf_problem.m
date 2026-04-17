function problem = prepare_scf_problem(sys, Eext, scfParams)
%PREPARE_SCF_PROBLEM Normalize SCF inputs for all solver backends.
%
% Returns a struct with consistent shapes, defaults, cached vectors, and
% active-space indexing for polarizable sites only.
%
% Fields:
%   nSites, polMask, alpha
%   activeSites, nPolSites, activeVecIdx
%   alpha_pol, alpha_pol_vec
%   Eext, Eext_vec, Eext_pol, Eext_pol_vec
%   mu0,  mu0_vec,  mu0_pol,  mu0_pol_vec
%   tol, maxIter, mixing, omega
%   verbose, printEvery, residualEvery, stopMetric

    if nargin < 3 || isempty(scfParams)
        scfParams = struct();
    end

    io.assert_atomic_units(sys);

    nSites = sys.n_sites;

    if ~isequal(size(Eext), [nSites, 3])
        error('Eext must be N x 3.');
    end

    polMask = logical(sys.site_is_polarizable(:));
    alpha   = sys.site_alpha(:);

    if numel(polMask) ~= nSites || numel(alpha) ~= nSites
        error('sys.site_is_polarizable and sys.site_alpha must both have length nSites.');
    end

    % ----------------------------
    % Solver defaults
    % ----------------------------
    tol = 1e-8;
    if isfield(scfParams, 'tol') && ~isempty(scfParams.tol)
        tol = scfParams.tol;
    end
    if ~isscalar(tol) || tol <= 0
        error('scfParams.tol must be a positive scalar.');
    end

    maxIter = 500;
    if isfield(scfParams, 'maxIter') && ~isempty(scfParams.maxIter)
        maxIter = scfParams.maxIter;
    end
    if ~isscalar(maxIter) || maxIter < 1
        error('scfParams.maxIter must be a positive integer-like scalar.');
    end

    mixing = 0.5;
    if isfield(scfParams, 'mixing') && ~isempty(scfParams.mixing)
        mixing = scfParams.mixing;
    end

    omega = 1.0;
    if isfield(scfParams, 'omega') && ~isempty(scfParams.omega)
        omega = scfParams.omega;
    end
    if omega <= 0 || omega >= 2
        error('scfParams.omega must satisfy 0 < omega < 2.');
    end

    verbose = false;
    if isfield(scfParams, 'verbose') && ~isempty(scfParams.verbose)
        verbose = logical(scfParams.verbose);
    end

    printEvery = 10;
    if isfield(scfParams, 'printEvery') && ~isempty(scfParams.printEvery)
        printEvery = scfParams.printEvery;
    end
    if ~isscalar(printEvery) || printEvery < 1
        error('scfParams.printEvery must be a positive scalar.');
    end

    residualEvery = 25;
    if isfield(scfParams, 'residualEvery') && ~isempty(scfParams.residualEvery)
        residualEvery = scfParams.residualEvery;
    end
    if ~isscalar(residualEvery) || residualEvery < 1
        error('scfParams.residualEvery must be a positive scalar.');
    end

    stopMetric = "max_dmu";
    if isfield(scfParams, 'stopMetric') && ~isempty(scfParams.stopMetric)
        stopMetric = lower(string(scfParams.stopMetric));
    end

    % ----------------------------
    % Initial induced dipoles
    % ----------------------------
    if isfield(scfParams, 'initial_mu') && ~isempty(scfParams.initial_mu)
        mu0 = scfParams.initial_mu;
        if ~isequal(size(mu0), [nSites, 3])
            error('scfParams.initial_mu must be N x 3.');
        end
    else
        mu0 = zeros(nSites, 3);
    end

    % ----------------------------
    % Active-space indexing
    % ----------------------------
    activeSites = find(polMask);
    nPolSites   = numel(activeSites);

    if nPolSites == 0
        error('prepare_scf_problem:NoPolarizableSites', ...
            'No polarizable sites are available for the SCF problem.');
    end

    activeVecIdx = zeros(3*nPolSites, 1);
    for k = 1:nPolSites
        i = activeSites(k);
        activeVecIdx(3*k-2:3*k) = util.block3(i);
    end

    alpha_pol = alpha(activeSites);

    Eext_pol = Eext(activeSites, :);
    mu0_pol  = mu0(activeSites, :);

    % xyz-stacked active-space vectors
    alpha_pol_vec = repelem(alpha_pol, 3, 1);
    Eext_pol_vec  = util.stack_xyz(Eext_pol);
    mu0_pol_vec   = util.stack_xyz(mu0_pol);

    % ----------------------------
    % Package result
    % ----------------------------
    problem = struct();

    % full-space data
    problem.nSites    = nSites;
    problem.polMask   = polMask;
    problem.alpha     = alpha;
    problem.Eext      = Eext;
    problem.Eext_vec  = util.stack_xyz(Eext);
    problem.mu0       = mu0;
    problem.mu0_vec   = util.stack_xyz(mu0);

    % active-space data
    problem.activeSites    = activeSites;
    problem.nPolSites      = nPolSites;
    problem.activeVecIdx   = activeVecIdx;

    problem.alpha_pol      = alpha_pol;
    problem.alpha_pol_vec  = alpha_pol_vec;

    problem.Eext_pol       = Eext_pol;
    problem.Eext_pol_vec   = Eext_pol_vec;

    problem.mu0_pol        = mu0_pol;
    problem.mu0_pol_vec    = mu0_pol_vec;

    % solver controls
    problem.tol            = tol;
    problem.maxIter        = maxIter;
    problem.mixing         = mixing;
    problem.omega          = omega;
    problem.verbose        = verbose;
    problem.printEvery     = printEvery;
    problem.residualEvery  = residualEvery;
    problem.stopMetric     = char(stopMetric);
end