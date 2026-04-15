function problem = prepare_scf_problem(sys, Eext, scfParams)
%PREPARE_SCF_PROBLEM Normalize SCF inputs for all solver backends.
%
% Returns a struct with consistent shapes, defaults, and cached vectors.
%
% Fields:
%   nSites, polMask, alpha
%   Eext, Eext_vec
%   mu0, mu0_vec
%   tol, maxIter, mixing, omega
%   verbose, printEvery
%   residualEvery
%   stopMetric

    if nargin < 3 || isempty(scfParams)
        scfParams = struct();
    end

    nSites = sys.n_sites;

    if ~isequal(size(Eext), [nSites, 3])
        error('Eext must be N x 3.');
    end

    polMask = logical(sys.site_is_polarizable(:));
    alpha   = sys.site_alpha(:);

    if numel(polMask) ~= nSites || numel(alpha) ~= nSites
        error('sys.site_is_polarizable and sys.site_alpha must both have length nSites.');
    end

    tol = 1e-8;
    if isfield(scfParams, 'tol') && ~isempty(scfParams.tol)
        tol = scfParams.tol;
    end

    maxIter = 500;
    if isfield(scfParams, 'maxIter') && ~isempty(scfParams.maxIter)
        maxIter = scfParams.maxIter;
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

    residualEvery = 25;
    if isfield(scfParams, 'residualEvery') && ~isempty(scfParams.residualEvery)
        residualEvery = scfParams.residualEvery;
    end

    stopMetric = 'max_dmu';
    if isfield(scfParams, 'stopMetric') && ~isempty(scfParams.stopMetric)
        stopMetric = lower(string(scfParams.stopMetric));
    end

    if isfield(scfParams, 'initial_mu') && ~isempty(scfParams.initial_mu)
        mu0 = scfParams.initial_mu;
        if ~isequal(size(mu0), [nSites, 3])
            error('scfParams.initial_mu must be N x 3.');
        end
    else
        mu0 = zeros(nSites, 3);
    end

    problem = struct();
    problem.nSites = nSites;
    problem.polMask = polMask;
    problem.alpha = alpha;

    problem.Eext = Eext;
    problem.Eext_vec = util.stack_xyz(Eext);

    problem.mu0 = mu0;
    problem.mu0_vec = util.stack_xyz(mu0);

    problem.tol = tol;
    problem.maxIter = maxIter;
    problem.mixing = mixing;
    problem.omega = omega;

    problem.verbose = verbose;
    problem.printEvery = printEvery;
    problem.residualEvery = residualEvery;
    problem.stopMetric = char(stopMetric);
end