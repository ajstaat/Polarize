function [mu, scf] = solve_scf_iterative(sys, Eext, scfParams)
%SOLVE_SCF_ITERATIVE Solve induced dipoles by fixed-point SCF iteration.
%
% Inputs
%   sys        struct with fields:
%              .site_is_polarizable   N x 1 logical
%              .site_alpha            N x 1 scalar polarizabilities
%
%   Eext       N x 3 external field
%
%   scfParams  struct with optional fields:
%              .tol            default 1e-8
%              .maxIter        default 500
%              .mixing         default 0.5
%              .initial_mu     optional N x 3
%              .softening      default 0
%              .use_thole      default true
%
% Outputs
%   mu         N x 3 induced dipoles
%   scf        struct with convergence metadata

if ~isfield(sys, 'site_is_polarizable') || isempty(sys.site_is_polarizable)
    error('sys.site_is_polarizable is missing or empty.');
end

if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
    error('sys.site_alpha is missing or empty.');
end

nSites = sys.n_sites;

if ~isequal(size(Eext), [nSites, 3])
    error('Eext must be N x 3.');
end

tol = 1e-8;
if isfield(scfParams, 'tol')
    tol = scfParams.tol;
end

maxIter = 500;
if isfield(scfParams, 'maxIter')
    maxIter = scfParams.maxIter;
end

mixing = 0.5;
if isfield(scfParams, 'mixing')
    mixing = scfParams.mixing;
end

softening = 0.0;
if isfield(scfParams, 'softening')
    softening = scfParams.softening;
end

use_thole = true;
if isfield(scfParams, 'use_thole')
    use_thole = scfParams.use_thole;
end

if isfield(scfParams, 'initial_mu') && ~isempty(scfParams.initial_mu)
    mu = scfParams.initial_mu;
    if ~isequal(size(mu), [nSites, 3])
        error('scfParams.initial_mu must be N x 3.');
    end
else
    mu = zeros(nSites, 3);
end

polMask = logical(sys.site_is_polarizable(:));
alpha = sys.site_alpha(:);

if numel(alpha) ~= nSites
    error('sys.site_alpha must have length N.');
end

history = zeros(maxIter, 1);
converged = false;

dipoleParams = struct();
dipoleParams.exclude_self = true;
dipoleParams.softening = softening;
dipoleParams.target_mask = polMask;
dipoleParams.source_mask = polMask;

for iter = 1:maxIter
    if use_thole
        Edip = thole.induced_field_from_dipoles_thole(sys, mu, dipoleParams);
    else
        Edip = thole.induced_field_from_dipoles(sys, mu, dipoleParams);
    end

    Etot = Eext + Edip;

    mu_new = zeros(nSites, 3);
    mu_new(polMask, :) = alpha(polMask) .* Etot(polMask, :);

    mu_mixed = (1 - mixing) * mu + mixing * mu_new;

    delta = mu_mixed - mu;
    err = max(sqrt(sum(delta.^2, 2)));
    history(iter) = err;

    mu = mu_mixed;

    if err < tol
        converged = true;
        history = history(1:iter);
        break;
    end
end

if ~converged
    history = history(1:maxIter);
end

scf = struct();
scf.converged = converged;
scf.nIter = numel(history);
scf.history = history;
scf.used_thole = use_thole;

end