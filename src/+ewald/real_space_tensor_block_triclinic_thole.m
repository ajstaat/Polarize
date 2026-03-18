function Tij = real_space_tensor_block_triclinic_thole(xvec, alpha, alpha_i, alpha_j, thole_a)
%REAL_SPACE_TENSOR_BLOCK_TRICLINIC_THOLE
% Real-space Ewald dipole tensor block with Thole correction.
%
% Inputs
%   xvec     1x3 Cartesian displacement vector
%   alpha    Ewald screening parameter
%   alpha_i  scalar polarizability of site i
%   alpha_j  scalar polarizability of site j
%   thole_a  scalar Thole parameter
%
% Output
%   Tij      3x3 damped real-space tensor block
%
% Form:
%   Tbare   = B I - C xx^T
%   dT      = (l3/x^3) I - (3 l5/x^5) xx^T
%   Tdamped = Tbare + dT

if numel(xvec) ~= 3
    error('xvec must be a 1x3 vector.');
end
if alpha <= 0
    error('alpha must be positive.');
end

Tbare = ewald.real_space_tensor_block_triclinic(xvec, alpha);
dT = ewald.real_space_tensor_block_triclinic_thole_correction(xvec, alpha_i, alpha_j, thole_a);

Tij = Tbare + dT;

end