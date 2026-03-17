function lattice = lattice_vectors_from_cell(cellpar)
%LATTICE_VECTORS_FROM_CELL Build 3x3 lattice matrix from cell parameters.
%
% Input
%   cellpar : [a b c alpha beta gamma]
%             lengths in any consistent length unit
%             angles in degrees
%
% Output
%   lattice : 3x3 matrix whose rows are lattice vectors:
%             [ax ay az
%              bx by bz
%              cx cy cz]
%
% Convention
%   a is along x
%   b lies in the xy plane
%   c has general triclinic orientation

if numel(cellpar) ~= 6
    error('cellpar must be a 1x6 vector: [a b c alpha beta gamma].');
end

a = cellpar(1);
b = cellpar(2);
c = cellpar(3);

alpha = deg2rad(cellpar(4));
beta  = deg2rad(cellpar(5));
gamma = deg2rad(cellpar(6));

% Basic checks
if a <= 0 || b <= 0 || c <= 0
    error('Cell lengths a, b, c must be positive.');
end

sinGamma = sin(gamma);
if abs(sinGamma) < 1e-14
    error('gamma is too close to 0 or 180 degrees; lattice is singular.');
end

% Standard triclinic construction
avec = [a, 0, 0];

bvec = [b*cos(gamma), ...
        b*sin(gamma), ...
        0];

cx = c*cos(beta);
cy = c*(cos(alpha) - cos(beta)*cos(gamma)) / sinGamma;

cz2 = c^2 - cx^2 - cy^2;
if cz2 < -1e-10
    error('Cell parameters are inconsistent; computed c_z^2 is negative.');
end
cz = sqrt(max(cz2, 0));

cvec = [cx, cy, cz];

lattice = [avec; bvec; cvec];

end