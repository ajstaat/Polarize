function tf = thole_f3f5_factors(r, alpha_i, alpha_j, a)
%THOLE_F3F5_FACTORS Compute Thole damping factors f3/f5 and corrections l3/l5.
%
% Inputs
%   r        scalar separation
%   alpha_i  scalar polarizability of site i
%   alpha_j  scalar polarizability of site j
%   a        scalar Thole parameter
%
% Output
%   tf       struct with fields:
%            .s
%            .f3
%            .f5
%            .l3
%            .l5
%
% Definitions
%   s  = a * ( r / (alpha_i*alpha_j)^(1/6) )^3
%   f3 = 1 - exp(-s)
%   f5 = 1 - (1 + s) exp(-s)
%   l3 = f3 - 1 = -exp(-s)
%   l5 = f5 - 1 = -(1 + s) exp(-s)

if r < 0
    error('r must be nonnegative.');
end
if alpha_i < 0 || alpha_j < 0
    error('Polarizabilities must be nonnegative.');
end
if a < 0
    error('Thole parameter a must be nonnegative.');
end

tf = struct();
tf.s = 0.0;
tf.f3 = 1.0;
tf.f5 = 1.0;
tf.l3 = 0.0;
tf.l5 = 0.0;

if r == 0 || alpha_i == 0 || alpha_j == 0 || a == 0
    return;
end

aij = (alpha_i * alpha_j)^(1/6);

if aij == 0
    return;
end

s = a * (r / aij)^3;
e = exp(-s);

tf.s = s;
tf.f3 = 1 - e;
tf.f5 = 1 - (1 + s) * e;
tf.l3 = -e;
tf.l5 = -(1 + s) * e;

end