function [Utotal, parts, meta] = dipole_energy_periodic_triclinic(mu, r, H, alpha, rcut, kcut, boundary)
%DIPOLE_ENERGY_PERIODIC_TRICLINIC
% Reference 3D Ewald summation for point dipoles in a triclinic cell.

if nargin < 7 || isempty(boundary)
    boundary = 'tinfoil';
end

if size(mu,2) ~= 3 || size(r,2) ~= 3
    error('mu and r must both be N x 3 arrays.');
end
if size(mu,1) ~= size(r,1)
    error('mu and r must have the same number of rows.');
end
if ~isequal(size(H), [3 3])
    error('H must be a 3x3 cell matrix with lattice vectors as columns.');
end
if alpha <= 0 || rcut <= 0 || kcut <= 0
    error('alpha, rcut, and kcut must be positive.');
end

boundary = lower(boundary);
if ~ismember(boundary, {'tinfoil','vacuum'})
    error('boundary must be ''tinfoil'' or ''vacuum''.');
end

V = abs(det(H));
if V <= 0
    error('Cell matrix must have positive nonzero volume.');
end

[Ureal, realMeta] = ewald.real_space_energy_triclinic(mu, r, H, alpha, rcut);
[Urecip, recipMeta] = ewald.reciprocal_space_energy_triclinic(mu, r, H, alpha, kcut);
Uself = ewald.self_term_dipole(mu, alpha);
Usurf = ewald.surface_term_dipole(mu, H, boundary);

Utotal = Ureal + Urecip + Uself + Usurf;

parts = struct();
parts.real = Ureal;
parts.recip = Urecip;
parts.self = Uself;
parts.surf = Usurf;

meta = struct();
meta.real_image_bounds = realMeta.real_image_bounds;
meta.num_kvec = recipMeta.num_kvec;
meta.volume = V;
meta.boundary = boundary;
meta.kcut = kcut;
meta.rcut = rcut;
meta.alpha = alpha;
end