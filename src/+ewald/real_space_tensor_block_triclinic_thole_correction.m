function dTij = real_space_tensor_block_triclinic_thole_correction(xvec, alpha_i, alpha_j, thole_a)
%REAL_SPACE_TENSOR_BLOCK_TRICLINIC_THOLE_CORRECTION
% Thole correction added to the bare Ewald real-space tensor.
%
% Inputs
%   xvec     1x3 Cartesian displacement vector
%   alpha_i  scalar polarizability of site i
%   alpha_j  scalar polarizability of site j
%   thole_a  scalar Thole parameter
%
% Output
%   dTij     3x3 correction such that
%            Tdamped = Tbare + dTij
%
% Correction form:
%   dTij = (l3/x^3) I - (3*l5/x^5) xx^T
%
% where
%   l3 = f3 - 1
%   l5 = f5 - 1

if numel(xvec) ~= 3
    error('xvec must be a 1x3 vector.');
end
if alpha_i < 0 || alpha_j < 0
    error('alpha_i and alpha_j must be nonnegative.');
end
if thole_a < 0
    error('thole_a must be nonnegative.');
end

x = norm(xvec);

if x == 0
    dTij = zeros(3,3);
    return;
end

tf = thole.thole_f3f5_factors(x, alpha_i, alpha_j, thole_a);

I3 = eye(3);
xxT = xvec(:) * xvec(:).';

dTij = (tf.l3 / x^3) * I3 - (3 * tf.l5 / x^5) * xxT;

end