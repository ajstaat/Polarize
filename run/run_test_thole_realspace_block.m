clear; clc;

xvec = [1.2, -0.7, 2.1];
alpha = 0.35;

alpha_i = 10.0;
alpha_j = 8.0;
thole_a = 0.39;

Tbare = ewald.real_space_tensor_block_triclinic(xvec, alpha);
Tdamp = ewald.real_space_tensor_block_triclinic_thole(xvec, alpha, alpha_i, alpha_j, thole_a);
dT = ewald.real_space_tensor_block_triclinic_thole_correction(xvec, alpha_i, alpha_j, thole_a);

fprintf('||Tdamp - (Tbare + dT)|| = %.16e\n', norm(Tdamp - (Tbare + dT), 'fro'));

tf = thole.thole_f3f5_factors(norm(xvec), alpha_i, alpha_j, thole_a);
disp(tf)