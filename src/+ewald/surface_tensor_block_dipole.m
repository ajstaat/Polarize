function Tsurf = surface_tensor_block_dipole(H, boundary)
%SURFACE_TENSOR_BLOCK_DIPOLE Surface tensor block for dipole Ewald operator.

V = abs(det(H));
if V <= 0
    error('Cell matrix must have positive nonzero volume.');
end

boundary = lower(boundary);

switch boundary
    case 'tinfoil'
        Tsurf = zeros(3,3);
    case 'vacuum'
        Tsurf = (4*pi / (3*V)) * eye(3);
    otherwise
        error('boundary must be ''tinfoil'' or ''vacuum''.');
end

end