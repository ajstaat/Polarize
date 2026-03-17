clear; clc;

mu = [0 1 0;
      0 1 0];

r = [0 0 0;
     0 0 10];

H0 = [ 1.00   0.35   0.20;
       0.00   0.95   0.25;
       0.00   0.00   1.10 ];

H0 = 40 * H0;

fprintf('Triclinic large-cell convergence test\n');
fprintf('-------------------------------------\n');

for s = [2 2.5 3 4 5]
    H = s * H0;

    p = ewald.choose_ewald_params(H, 1e-8, 'tinfoil', 0.9);
    [U, parts, meta] = ewald.dipole_energy_periodic_triclinic( ...
        mu, r, H, p.alpha, p.rcut, p.kcut, p.boundary);

    fprintf(['scale = %4.1f   V = %12.4f   U = %+ .12e   ', ...
             'err = %.3e   nk = %d\n'], ...
             s, meta.volume, U, abs(U - bare_pair_energy(mu, r)), meta.num_kvec);
end

Ubare = bare_pair_energy(mu, r);
fprintf('\nBare pair energy = %+ .12e\n', Ubare);

fprintf('\nOrthorhombic special-case check (using triclinic code)\n');
fprintf('-----------------------------------------------------\n');

for L = [80 100 120 160 200]
    H = diag([L L L]);
    p = ewald.choose_ewald_params(H, 1e-8, 'tinfoil', 0.9);
    [U, parts, meta] = ewald.dipole_energy_periodic_triclinic( ...
        mu, r, H, p.alpha, p.rcut, p.kcut, p.boundary);

    fprintf('L = %3d   U = %+ .12e   err = %.3e   nk = %d\n', ...
            L, U, abs(U - Ubare), meta.num_kvec);
end

fprintf('\nRotation invariance test\n');
fprintf('------------------------\n');

H = 3.0 * H0;
p = ewald.choose_ewald_params(H, 1e-8, 'tinfoil', 0.9);
[Uref, parts_ref, meta_ref] = ewald.dipole_energy_periodic_triclinic( ...
    mu, r, H, p.alpha, p.rcut, p.kcut, p.boundary);

fprintf('Reference energy = %+ .15e\n', Uref);

angles_deg = [20, 35, 50];

for t = 1:numel(angles_deg)
    ang = angles_deg(t) * pi / 180;

    switch t
        case 1
            axis = [1; 0; 0];
        case 2
            axis = [0; 1; 1];
        case 3
            axis = [1; 2; 3];
    end
    axis = axis / norm(axis);

    R = axis_angle_rotation(axis, ang);

    Hrot  = R * H;
    rrot  = (R * r.').';
    murot = (R * mu.').';

    prot = ewald.choose_ewald_params(Hrot, 1e-8, 'tinfoil', 0.9);
    [Urot, parts_rot, meta_rot] = ewald.dipole_energy_periodic_triclinic( ...
        murot, rrot, Hrot, prot.alpha, prot.rcut, prot.kcut, prot.boundary);

    fprintf('angle = %5.1f deg   Urot = %+ .15e   |dU| = %.3e\n', ...
            angles_deg(t), Urot, abs(Urot - Uref));
end

function Ubare = bare_pair_energy(mu, r)
    rij = r(1,:) - r(2,:);
    dist = norm(rij);

    Ubare = dot(mu(1,:), mu(2,:)) / dist^3 ...
          - 3 * dot(mu(1,:), rij) * dot(mu(2,:), rij) / dist^5;
end

function R = axis_angle_rotation(axis, theta)
    ax = axis(:) / norm(axis);

    x = ax(1);
    y = ax(2);
    z = ax(3);

    K = [  0  -z   y;
           z   0  -x;
          -y   x   0 ];

    R = eye(3) + sin(theta)*K + (1-cos(theta))*(K*K);
end