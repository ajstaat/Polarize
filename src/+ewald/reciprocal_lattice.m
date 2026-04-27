function G = reciprocal_lattice(H)
%RECIPROCAL_LATTICE Reciprocal lattice matrix for a triclinic cell.
%
% G = ewald.reciprocal_lattice(H)
%
% Input
%   H   3x3 direct lattice matrix, columns are lattice vectors
%
% Output
%   G   3x3 reciprocal lattice matrix, columns are reciprocal vectors
%
% Notes
%   - Uses the convention G = 2*pi*(H^{-1})^T, so that
%       a_i · g_j = 2*pi*delta_ij
%     when H columns are the direct lattice vectors.

    validateattributes(H, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
        mfilename, 'H', 1);

    V = det(H);
    if abs(V) < 1e-14
        error('ewald:reciprocal_lattice:SingularCell', ...
            'H is singular or nearly singular.');
    end

    G = 2 * pi * inv(H).';
end