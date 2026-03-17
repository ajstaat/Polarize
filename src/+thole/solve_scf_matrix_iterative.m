function [mu, scf] = solve_scf_matrix_iterative(sys, Eext, T, scfParams)
%SOLVE_SCF_MATRIX_ITERATIVE Solve SCF using explicit interaction matrix T.
%
% Inputs
%   sys        with .site_is_polarizable and .site_alpha
%   Eext       N x 3 external field
%   T          3N x 3N dipole interaction matrix
%   scfParams  same idea as solve_scf_iterative
%
% Outputs
%   mu         N x 3 induced dipoles
%   scf        convergence info

nSites = sys.n_sites;

if ~isequal(size(Eext), [nSites, 3])
    error('Eext must be N x 3.');
end

if ~isequal(size(T), [3*nSites, 3*nSites])
    error('T must be 3N x 3N.');
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

Eext_vec = util.stack_xyz(Eext);
mu_vec = util.stack_xyz(mu);

history = zeros(maxIter, 1);
converged = false;

for iter = 1:maxIter
    Edip_vec = T * mu_vec;
    Etot_vec = Eext_vec + Edip_vec;
    Etot = util.unstack_xyz(Etot_vec);

    mu_new = zeros(nSites, 3);
    mu_new(polMask, :) = alpha(polMask) .* Etot(polMask, :);

    mu_new_vec = util.stack_xyz(mu_new);
    mu_mixed_vec = (1 - mixing) * mu_vec + mixing * mu_new_vec;

    delta_vec = mu_mixed_vec - mu_vec;
    delta = util.unstack_xyz(delta_vec);
    err = max(sqrt(sum(delta.^2, 2)));
    history(iter) = err;

    mu_vec = mu_mixed_vec;

    if err < tol
        converged = true;
        history = history(1:iter);
        break;
    end
end

if ~converged
    history = history(1:maxIter);
end

mu = util.unstack_xyz(mu_vec);

scf = struct();
scf.converged = converged;
scf.nIter = numel(history);
scf.history = history;
scf.used_matrix_solver = true;

end