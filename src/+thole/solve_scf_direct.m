function [mu, direct] = solve_scf_direct(sys, Eext, T)
%SOLVE_SCF_DIRECT Solve induced dipoles by direct linear solve.
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
    Eext_vec = util.stack_xyz(Eext);

    M = eye(3*nSites) - A * T;
    b = A * Eext_vec;

    mu_vec = M \ b;
    mu = util.unstack_xyz(mu_vec);

    residual = M * mu_vec - b;

    direct = struct();
    direct.method = 'direct';
    direct.residual_norm = norm(residual);
    direct.used_matrix_solver = true;
end