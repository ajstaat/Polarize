function Usurf = surface_term_dipole(mu, H, boundary)
%SURFACE_TERM_DIPOLE Surface term for dipole Ewald sum.

boundary = lower(boundary);
V = abs(det(H));
M = sum(mu, 1);

switch boundary
    case 'tinfoil'
        Usurf = 0.0;
    case 'vacuum'
        Usurf = (2*pi / (3*V)) * dot(M, M);
    otherwise
        error('boundary must be ''tinfoil'' or ''vacuum''.');
end

end