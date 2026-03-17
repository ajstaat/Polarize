function Tself = self_tensor_block_dipole(alpha)
%SELF_TENSOR_BLOCK_DIPOLE Self tensor block for dipole Ewald operator.

if alpha <= 0
    error('alpha must be positive.');
end

Tself = -(4 * alpha^3 / (3 * sqrt(pi))) * eye(3);

end