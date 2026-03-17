function Tij = real_space_tensor_block_triclinic(xvec, alpha)
%REAL_SPACE_TENSOR_BLOCK_TRICLINIC Real-space Ewald dipole tensor block.
%
% Input
%   xvec   1x3 Cartesian displacement vector
%   alpha  Ewald screening parameter
%
% Output
%   Tij    3x3 tensor block

if numel(xvec) ~= 3
    error('xvec must be a 1x3 vector.');
end
if alpha <= 0
    error('alpha must be positive.');
end

x = norm(xvec);

if x == 0
    Tij = zeros(3,3);
    return;
end

erfcax = erfc(alpha * x);
expax2 = exp(-(alpha^2) * (x^2));

B = erfcax / x^3 + (2*alpha/sqrt(pi)) * expax2 / x^2;
C = 3*erfcax / x^5 + (2*alpha/sqrt(pi)) * ...
    (2*alpha^2 / x^2 + 3 / x^4) * expax2;

I3 = eye(3);
xxT = xvec(:) * xvec(:).';

Tij = B * I3 - C * xxT;
end