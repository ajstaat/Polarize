function [Tpol, opinfo] = assemble_nonperiodic_interaction_matrix(sys, problem, ewaldParams, scfParams)
%ASSEMBLE_NONPERIODIC_INTERACTION_MATRIX
% Build polarizable-only dipole interaction matrix for nonperiodic systems.
%
% Inputs
%   sys         struct with fields:
%                 .site_pos
%                 .site_alpha
%                 .site_is_polarizable
%                 .thole_a
%   problem     struct from thole.prepare_scf_problem(...)
%   ewaldParams unused here for now, included for API compatibility
%   scfParams   optional struct with fields:
%                 .softening   scalar, default 0
%                 .use_thole   logical, default true
%                 .rcut        scalar cutoff in bohr, optional
%                 .verbose     logical, default false
%                 .printEvery  positive integer, default 25
%
% Output
%   Tpol        3Np x 3Np interaction matrix such that:
%                 E_dip_pol_vec = Tpol * mu_pol_vec
%   opinfo      metadata struct describing active-space mapping
%
% Notes
%   - Only polarizable-to-polarizable couplings are included.
%   - Self blocks are zero.
%   - Matrix is assembled in polarizable-site space only.
%   - If rcut is supplied, only active-site pairs with bare separation
%     |r_ij| <= rcut are included.

    io.assert_atomic_units(sys);

    if nargin < 3
        ewaldParams = struct(); %#ok<NASGU>
    end
    if nargin < 4 || isempty(scfParams)
        scfParams = struct();
    end

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('sys.site_pos is missing or empty.');
    end
    if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
        error('sys.site_alpha is missing or empty.');
    end
    if ~isfield(sys, 'site_is_polarizable') || isempty(sys.site_is_polarizable)
        error('sys.site_is_polarizable is missing or empty.');
    end
    if ~isfield(sys, 'thole_a') || isempty(sys.thole_a)
        error('sys.thole_a is missing or empty.');
    end

    if nargin < 2 || isempty(problem)
        error('problem struct from thole.prepare_scf_problem is required.');
    end
    if ~isfield(problem, 'activeSites') || ~isfield(problem, 'nPolSites') ...
            || ~isfield(problem, 'activeVecIdx')
        error('problem must include activeSites, nPolSites, and activeVecIdx.');
    end

    nSites = sys.n_sites;
    sites  = problem.activeSites(:);
    nPol   = problem.nPolSites;

    if numel(sites) ~= nPol
        error('problem.activeSites and problem.nPolSites are inconsistent.');
    end

    opts = struct();
    opts.softening = 0.0;
    if isfield(scfParams, 'softening') && ~isempty(scfParams.softening)
        opts.softening = scfParams.softening;
    end

    opts.use_thole = true;
    if isfield(scfParams, 'use_thole') && ~isempty(scfParams.use_thole)
        opts.use_thole = scfParams.use_thole;
    end

    rcut2 = inf;
    if isfield(scfParams, 'rcut') && ~isempty(scfParams.rcut)
        rcut = scfParams.rcut;
        if ~isscalar(rcut) || rcut <= 0
            error('scfParams.rcut must be a positive scalar when provided.');
        end
        rcut2 = rcut^2;
    end

    verbose = false;
    if isfield(scfParams, 'verbose') && ~isempty(scfParams.verbose)
        verbose = logical(scfParams.verbose);
    end

    printEvery = 25;
    if isfield(scfParams, 'printEvery') && ~isempty(scfParams.printEvery)
        printEvery = scfParams.printEvery;
    end

    Tpol = zeros(3*nPol, 3*nPol);

    pos_pol   = sys.site_pos(sites, :);
    alpha_pol = sys.site_alpha(sites);

    if verbose
        fprintf('Nonperiodic interaction matrix: %d total sites, %d polarizable sites\n', ...
            nSites, nPol);
        fprintf('Compressed matrix size: %d x %d\n', 3*nPol, 3*nPol);
        if isfinite(rcut2)
            fprintf('Applying nonperiodic cutoff: rcut = %.6f bohr\n', sqrt(rcut2));
        end
    end

    tStart = tic;
    nPairsKept = 0;

    for a = 1:nPol
        Ia = (3*a-2):(3*a);
        ra = pos_pol(a, :);
        aa = alpha_pol(a);

        for b = (a+1):nPol
            rb = pos_pol(b, :);
            r2_bare = sum((ra - rb).^2);

            if r2_bare > rcut2
                continue;
            end

            Ib = (3*b-2):(3*b);

            Tij = thole.dipole_tensor_block( ...
                ra, ...
                rb, ...
                aa, alpha_pol(b), ...
                sys.thole_a, opts);

            Tpol(Ia, Ib) = Tij;
            Tpol(Ib, Ia) = Tij.';
            nPairsKept = nPairsKept + 1;
        end

        if verbose && (a == 1 || a == nPol || mod(a, printEvery) == 0)
            frac      = a / max(nPol, 1);
            elapsed   = toc(tStart);
            totalEst  = elapsed / max(frac, eps);
            remainEst = max(totalEst - elapsed, 0);

            fprintf([' assembling Tpol: active row %4d / %4d ' ...
                     '(site %4d / %4d, %.1f%%) elapsed %.1f s ETA %.1f s\n'], ...
                a, nPol, sites(a), nSites, 100*frac, elapsed, remainEst);
        end
    end

    if verbose
        fprintf('Finished assembling compressed nonperiodic interaction matrix in %.2f s\n', ...
            toc(tStart));
    end

    opinfo = struct();
    opinfo.space        = 'polarizable_only';
    opinfo.nSites       = nSites;
    opinfo.nPolSites    = nPol;
    opinfo.activeSites  = sites;
    opinfo.activeVecIdx = problem.activeVecIdx;
    opinfo.rcut         = sqrt(rcut2);
    opinfo.nPairBlocksKept = nPairsKept;
end