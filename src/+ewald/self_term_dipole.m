function Uself = self_term_dipole(mu, alpha)
%SELF_TERM_DIPOLE Dipole Ewald self term.
%
% Uself = ewald.self_term_dipole(mu, alpha)
%
% Inputs
%   mu     N x 3 dipole array
%   alpha  Ewald splitting parameter
%
% Output
%   Uself  scalar self-energy contribution
%
% Notes
%   - Uses the standard dipole Ewald self term:
%       Uself = -(2*alpha^3 / (3*sqrt(pi))) * sum_i |mu_i|^2

    validateattributes(mu, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
        mfilename, 'mu', 1);
    validateattributes(alpha, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'alpha', 2);

    mu2sum = sum(sum(mu.^2, 2));
    Uself = -(2 * alpha^3 / (3 * sqrt(pi))) * mu2sum;
end