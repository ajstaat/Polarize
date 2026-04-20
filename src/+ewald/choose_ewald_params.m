function params = choose_ewald_params(box, tol, boundary, rcut_fraction)
%CHOOSE_EWALD_PARAMS Choose practical triclinic Ewald parameters.
%
% params = ewald.choose_ewald_params(box)
% params = ewald.choose_ewald_params(box, tol)
% params = ewald.choose_ewald_params(box, tol, boundary)
% params = ewald.choose_ewald_params(box, tol, boundary, rcut_fraction)
%
% Inputs
%   box            3x3 direct lattice matrix, columns are lattice vectors
%   tol            target Ewald truncation tolerance, default 1e-6
%   boundary       boundary condition label, default 'tinfoil'
%   rcut_fraction  fraction of half the shortest box-vector norm to use for
%                  real-space cutoff, default 0.9
%
% Output
%   params         struct with fields:
%       .rcut
%       .alpha
%       .kcut
%       .tol
%       .boundary
%       .rcut_fraction
%       .Lmin
%
% Notes
%   - This is a practical heuristic, not an optimality guarantee.
%   - The choice mirrors the old main-branch periodic workflow while fitting
%     the refactored package structure.

    if nargin < 2 || isempty(tol)
        tol = 1e-6;
    end
    if nargin < 3 || isempty(boundary)
        boundary = 'tinfoil';
    end
    if nargin < 4 || isempty(rcut_fraction)
        rcut_fraction = 0.9;
    end

    validateattributes(box, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
        mfilename, 'box', 1);

    validateattributes(tol, {'double'}, {'scalar', 'real', 'finite', '>', 0, '<', 1}, ...
        mfilename, 'tol', 2);

    if ~(ischar(boundary) || isstring(boundary))
        error('ewald:choose_ewald_params:BadBoundary', ...
            'boundary must be a character vector or string scalar.');
    end
    boundary = char(string(boundary));

    validateattributes(rcut_fraction, {'double'}, ...
        {'scalar', 'real', 'finite', '>', 0, '<=', 1}, ...
        mfilename, 'rcut_fraction', 4);

    a = box(:,1);
    b = box(:,2);
    c = box(:,3);

    Lmin = min([norm(a), norm(b), norm(c)]);
    rcut = rcut_fraction * (Lmin / 2);

    if rcut <= 0
        error('ewald:choose_ewald_params:NonPositiveRcut', ...
            'Computed rcut must be positive.');
    end

    alpha = erfcinv(tol) / rcut;
    kcut  = 2 * alpha * sqrt(log(1 / tol));

    params = struct();
    params.rcut = rcut;
    params.alpha = alpha;
    params.kcut = kcut;
    params.tol = tol;
    params.boundary = boundary;
    params.rcut_fraction = rcut_fraction;
    params.Lmin = Lmin;
end