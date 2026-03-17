function [mu, direct] = solve_scf_direct(sys, Eext, T)
%SOLVE_SCF_DIRECT Solve induced dipoles by direct linear solve.
%
% Inputs
%   sys    working system struct
%   Eext   N x 3 external field
%   T      3N x 3N dipole interaction matrix
%
% Outputs
%   mu     N x 3 induced dipoles
%   direct struct with metadata
%
% Solves:
%   (I - A*T) mu = A * Eext
%
% where A is the block-diagonal polarizability matrix.

nSites = sys.n_sites;

if ~isequal(size(Eext), [nSites, 3])
    error('Eext must be N x 3.');
end

if ~isequal(size(T), [3*nSites, 3*nSites])
    error('T must be 3N x 3N.');
end

A = thole.build_alpha_matrix(sys);
I = eye(3*nSites);

Eext_vec = util.stack_xyz(Eext);

M = I - A*T;
b = A * Eext_vec;

mu_vec = M \ b;
mu = util.unstack_xyz(mu_vec);

residual = M * mu_vec - b;

direct = struct();
direct.method = 'direct';
direct.A = A;
direct.M = M;
direct.b = b;
direct.residual_norm = norm(residual);
direct.condition_estimate = cond(M);

end