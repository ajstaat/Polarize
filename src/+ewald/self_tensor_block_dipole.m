function Tself = self_tensor_block_dipole(alpha)
%SELF_TENSOR_BLOCK_DIPOLE Analytic Ewald self block for dipole FIELD operator.
%
% Tself = ewald.self_tensor_block_dipole(alpha)
%
% Input
%   alpha   Ewald screening parameter
%
% Output
%   Tself   3x3 self block in the convention used by the dipole field operator
%           assembled into Tpol for the SCF equation
%
% Notes
%   The Ewald dipole self energy is
%       Uself = -(2*alpha^3/(3*sqrt(pi))) * sum_i |mu_i|^2
%
%   In this code path, Tpol is used as the dipole FIELD operator in
%       (I - A*Tpol) mu = A*Eext
%
%   so the self block included in Tpol has the opposite sign from the
%   energy-tensor convention:
%       Tself = +(4*alpha^3/(3*sqrt(pi))) * I

    validateattributes(alpha, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'alpha', 1);

    coeff = (4 * alpha^3 / (3 * sqrt(pi)));
    Tself = coeff * eye(3);
end