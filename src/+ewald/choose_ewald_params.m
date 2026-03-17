function params = choose_ewald_params(box, tol, boundary, rcut_fraction)
%CHOOSE_EWALD_PARAMS Choose practical triclinic Ewald parameters.

if nargin < 2 || isempty(tol)
    tol = 1e-6;
end
if nargin < 3 || isempty(boundary)
    boundary = 'tinfoil';
end
if nargin < 4 || isempty(rcut_fraction)
    rcut_fraction = 0.9;
end

if tol <= 0 || tol >= 1
    error('tol must satisfy 0 < tol < 1.');
end
if rcut_fraction <= 0 || rcut_fraction > 1
    error('rcut_fraction must satisfy 0 < rcut_fraction <= 1.');
end
if ~isequal(size(box), [3 3])
    error('box must be a 3x3 cell matrix.');
end

a = box(:,1);
b = box(:,2);
c = box(:,3);

Lmin = min([norm(a), norm(b), norm(c)]);
rcut = rcut_fraction * (Lmin / 2);

alpha = erfcinv(tol) / rcut;
kcut = 2 * alpha * sqrt(log(1/tol));

params = struct();
params.rcut = rcut;
params.alpha = alpha;
params.kcut = kcut;
params.tol = tol;
params.boundary = boundary;
params.rcut_fraction = rcut_fraction;
end