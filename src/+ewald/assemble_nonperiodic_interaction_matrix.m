function T = assemble_nonperiodic_interaction_matrix(sys, ewaldParams, scfParams)
%ASSEMBLE_NONPERIODIC_INTERACTION_MATRIX Build 3N x 3N dipole interaction matrix.
%
% Inputs
%   sys         struct with fields:
%               .site_pos
%               .site_alpha
%               .site_is_polarizable
%               .thole_a
%
%   ewaldParams unused here for now, included for API compatibility
%
%   scfParams   optional struct with fields:
%               .softening   scalar, default 0
%               .use_thole   logical, default true
%               .verbose     logical, default false
%               .printEvery  positive integer, default 25
%
% Output
%   T           3N x 3N interaction matrix such that:
%               E_dip_vec = T * mu_vec
%
% Notes
%   - Only polarizable-to-polarizable couplings are included.
%   - Self blocks are zero.

if nargin < 2
    ewaldParams = struct(); %#ok<NASGU>
end

if nargin < 3 || isempty(scfParams)
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

nSites = sys.n_sites;
T = zeros(3*nSites, 3*nSites);

polMask = logical(sys.site_is_polarizable(:));
alpha = sys.site_alpha(:);

opts = struct();
opts.softening = 0.0;
if isfield(scfParams, 'softening')
    opts.softening = scfParams.softening;
end

opts.use_thole = true;
if isfield(scfParams, 'use_thole')
    opts.use_thole = scfParams.use_thole;
end

verbose = false;
if isfield(scfParams, 'verbose') && ~isempty(scfParams.verbose)
    verbose = logical(scfParams.verbose);
end

printEvery = 25;
if isfield(scfParams, 'printEvery') && ~isempty(scfParams.printEvery)
    printEvery = scfParams.printEvery;
end

nPol = nnz(polMask);
if verbose
    fprintf('Nonperiodic interaction matrix: %d total sites, %d polarizable sites\n', ...
        nSites, nPol);
    fprintf('Matrix size: %d x %d\n', 3*nSites, 3*nSites);
end

tStart = tic;
nPolDone = 0;

for i = 1:nSites
    if ~polMask(i)
        continue;
    end

    nPolDone = nPolDone + 1;
    ii = util.block3(i);

    for j = 1:nSites
        if ~polMask(j)
            continue;
        end

        if i == j
            continue;
        end

        jj = util.block3(j);

        Tij = thole.dipole_tensor_block( ...
            sys.site_pos(i, :), ...
            sys.site_pos(j, :), ...
            alpha(i), alpha(j), ...
            sys.thole_a, opts);

        T(ii, jj) = Tij;
    end

    if verbose && (nPolDone == 1 || nPolDone == nPol || mod(nPolDone, printEvery) == 0)
        frac = nPolDone / max(nPol, 1);
        elapsed = toc(tStart);
        totalEst = elapsed / max(frac, eps);
        remainEst = max(totalEst - elapsed, 0);

        fprintf(['  assembling T: polarizable row %4d / %4d  ' ...
                 '(site %4d / %4d, %.1f%%)  elapsed %.1f s  ETA %.1f s\n'], ...
                 nPolDone, nPol, i, nSites, 100*frac, elapsed, remainEst);
    end
end

if verbose
    fprintf('Finished assembling nonperiodic interaction matrix in %.2f s\n', toc(tStart));
end

end