function frame = compute_molecule_frame(X, varargin)
%COMPUTE_MOLECULE_FRAME Build a canonical local frame for a mostly planar molecule.
%
% frame = builder.compute_molecule_frame(X)
% frame = builder.compute_molecule_frame(X, ...)
%
% Inputs
%   X   N x 3 Cartesian coordinates for one molecule
%
% Optional name-value inputs
%   'ReferenceAxis'   1x3 vector used to choose the sign of the plane normal.
%                     Default = [0 1 0] (positive b direction)
%   'Verbose'         logical, default false
%
% Output
%   frame struct with fields:
%       .com
%       .e1
%       .e2
%       .e3
%       .local_coords
%       .farthest_atom_index
%
% Convention
%   - e3 is the best-fit plane normal, with sign chosen so that
%       dot(e3, ReferenceAxis) >= 0
%   - e1 lies in the molecular plane and points from COM toward the farthest atom
%     after projection into the plane
%   - e2 = cross(e3, e1)

    p = inputParser;
    addRequired(p, 'X', @(x) isnumeric(x) && size(x,2) == 3);
    addParameter(p, 'ReferenceAxis', [0 1 0], @(x) isnumeric(x) && numel(x) == 3);
    addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
    parse(p, X, varargin{:});
    opt = p.Results;

    if size(X,1) < 3
        error('builder:compute_molecule_frame:TooFewAtoms', ...
            'Need at least 3 atoms to define a molecular frame.');
    end

    com = mean(X, 1);
    Xc = X - com;

    % Best-fit plane normal from SVD
    [~, ~, V] = svd(Xc, 'econ');
    e3 = V(:, end).';

    refAxis = reshape(opt.ReferenceAxis, 1, 3);
    nref = norm(refAxis);
    if nref == 0
        error('builder:compute_molecule_frame:BadReferenceAxis', ...
            'ReferenceAxis must have nonzero norm.');
    end
    refAxis = refAxis / nref;

    % Fix sign of plane normal
    if dot(e3, refAxis) < 0
        e3 = -e3;
    end

    % Choose in-plane axis from farthest atom from COM
    r2 = sum(Xc.^2, 2);
    [~, iFar] = max(r2);
    v = Xc(iFar, :);

    % Project into molecular plane
    v = v - dot(v, e3) * e3;
    nv = norm(v);

    if nv < 1e-12
        error('builder:compute_molecule_frame:DegenerateInPlaneAxis', ...
            'Failed to construct an in-plane axis from the farthest atom.');
    end

    e1 = v / nv;
    e2 = cross(e3, e1);
    e2 = e2 / norm(e2);

    % Re-orthogonalize e1 for numerical cleanliness
    e1 = cross(e2, e3);
    e1 = e1 / norm(e1);

    % Local coordinates in canonical frame
    local_x = Xc * e1.';
    local_y = Xc * e2.';
    local_z = Xc * e3.';
    local_coords = [local_x, local_y, local_z];

    frame = struct();
    frame.com = com;
    frame.e1 = e1;
    frame.e2 = e2;
    frame.e3 = e3;
    frame.local_coords = local_coords;
    frame.farthest_atom_index = iFar;

    if opt.Verbose
        fprintf('Molecule frame summary:\n');
        fprintf('  COM                    = [%9.4f %9.4f %9.4f]\n', com);
        fprintf('  e1                     = [%9.4f %9.4f %9.4f]\n', e1);
        fprintf('  e2                     = [%9.4f %9.4f %9.4f]\n', e2);
        fprintf('  e3 (normal)            = [%9.4f %9.4f %9.4f]\n', e3);
        fprintf('  farthest atom index    = %d\n', iFar);
    end
end