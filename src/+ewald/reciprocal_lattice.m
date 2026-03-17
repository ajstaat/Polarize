function G = reciprocal_lattice(H)
%RECIPROCAL_LATTICE Reciprocal lattice matrix for a triclinic cell.
%
% Input
%   H   3x3 direct lattice matrix, columns are lattice vectors
%
% Output
%   G   3x3 reciprocal lattice matrix, columns are reciprocal vectors

if ~isequal(size(H), [3 3])
    error('H must be 3x3.');
end

V = det(H);
if abs(V) < 1e-14
    error('H is singular or nearly singular.');
end

G = 2*pi * inv(H).';

end