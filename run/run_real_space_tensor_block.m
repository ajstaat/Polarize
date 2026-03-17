clear; clc;

xvec = [1.2, -0.7, 2.1];
alpha = 0.35;

mui = [0.3, -0.4, 0.8];
muj = [-0.2, 0.5, 0.1];

Tij = ewald.real_space_tensor_block_triclinic(xvec, alpha);

% Tensor contraction
Utensor = dot(mui, Tij * muj.');

% Scalar reference formula
x = norm(xvec);
erfcax = erfc(alpha * x);
expax2 = exp(-(alpha^2) * (x^2));

B = erfcax / x^3 + (2*alpha/sqrt(pi)) * expax2 / x^2;
C = 3*erfcax / x^5 + (2*alpha/sqrt(pi)) * ...
    (2*alpha^2 / x^2 + 3 / x^4) * expax2;

Uscalar = B * dot(mui, muj) - C * dot(mui, xvec) * dot(muj, xvec);

fprintf('Utensor = %+ .15e\n', Utensor);
fprintf('Uscalar = %+ .15e\n', Uscalar);
fprintf('|diff|   = %.3e\n', abs(Utensor - Uscalar));