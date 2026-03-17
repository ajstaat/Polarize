function f = thole_scale_factor(r, alpha_i, alpha_j, a)
%THOLE_SCALE_FACTOR Scalar Thole damping factor for dipole coupling.
%
% Inputs
%   r        scalar separation
%   alpha_i  scalar polarizability of site i
%   alpha_j  scalar polarizability of site j
%   a        scalar Thole parameter
%
% Output
%   f        scalar damping factor in [0,1]
%
% Uses a simple exponential damping form:
%   u = r / (alpha_i * alpha_j)^(1/6)
%   f = 1 - exp(-a * u^3)

if alpha_i < 0 || alpha_j < 0
    error('Polarizabilities must be nonnegative.');
end

if r < 0
    error('Distance r must be nonnegative.');
end

if alpha_i == 0 || alpha_j == 0
    f = 1.0;
    return;
end

alpha_ij = (alpha_i * alpha_j)^(1/6);

if alpha_ij == 0
    f = 1.0;
    return;
end

u = r / alpha_ij;
f = 1 - exp(-a * u^3);

% Numerical safety
f = max(0.0, min(1.0, f));

end