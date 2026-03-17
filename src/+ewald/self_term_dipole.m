function Uself = self_term_dipole(mu, alpha)
%SELF_TERM_DIPOLE Dipole Ewald self term.

mu2sum = sum(sum(mu.^2, 2));
Uself = -(2 * alpha^3 / (3 * sqrt(pi))) * mu2sum;

end