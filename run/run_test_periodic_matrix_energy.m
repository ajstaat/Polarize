clear; clc;

mu = [0 1 0;
      0 1 0];

r = [0 0 0;
     0 0 10];

H0 = [ 1.00   0.35   0.20;
       0.00   0.95   0.25;
       0.00   0.00   1.10 ];
H = 3.0 * 40 * H0;

p = ewald.choose_ewald_params(H, 1e-8, 'tinfoil', 0.9);

% Trusted energy path
[Uref, parts_ref, meta_ref] = ewald.dipole_energy_periodic_triclinic( ...
    mu, r, H, p.alpha, p.rcut, p.kcut, p.boundary);

% Build minimal sys for matrix path
sys = struct();
sys.site_pos = r;
sys.n_sites = size(r,1);
sys.site_is_polarizable = true(sys.n_sites,1);

% Note: assemble_periodic_interaction_matrix expects super_lattice or lattice
% with row-vectors, but triclinic Ewald code uses column-vectors.
% So pass row-vector version here.
sys.super_lattice = H.';

[T, parts_mat, meta_mat] = ewald.assemble_periodic_interaction_matrix(sys, p, struct());

mu_vec = util.stack_xyz(mu);
Umat = 0.5 * (mu_vec.' * T * mu_vec);

fprintf('Reference energy = %+ .15e\n', Uref);
fprintf('Matrix energy    = %+ .15e\n', Umat);
fprintf('|diff|           = %.3e\n', abs(Uref - Umat));