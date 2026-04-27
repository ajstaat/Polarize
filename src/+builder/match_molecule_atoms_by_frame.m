function match = match_molecule_atoms_by_frame(X_template, species_template, X_target, species_target, varargin)
%MATCH_MOLECULE_ATOMS_BY_FRAME Match atoms between two corresponding molecules
% using canonical local molecular frames.
%
% match = builder.match_molecule_atoms_by_frame(X_template, species_template, X_target, species_target)
% match = builder.match_molecule_atoms_by_frame(..., 'ReferenceAxis', bvec, ...)
%
% Inputs
%   X_template        N x 3 Cartesian coordinates for template molecule
%   species_template  N x 1 cell array of element symbols for template
%   X_target          N x 3 Cartesian coordinates for target molecule
%   species_target    N x 1 cell array of element symbols for target
%
% Optional name-value inputs
%   'ReferenceAxis'   1x3 vector used to fix the sign of the plane normal.
%                     Default = [0 1 0]
%   'DistanceTol'     scalar local-coordinate matching tolerance, default 1e-3
%   'Verbose'         logical, default false
%
% Output
%   match struct with fields:
%       .template_frame
%       .target_frame
%       .template_to_target      N x 1 integer map
%       .target_to_template      N x 1 integer map
%       .template_local_coords   N x 3
%       .target_local_coords     N x 3
%       .pair_distance           N x 1 matched local-frame distances
%       .max_pair_distance
%       .rms_pair_distance
%       .is_within_tolerance
%
% Notes
%   - Both molecules must contain the same multiset of elements.
%   - Matching is performed element-by-element.
%   - This is intended for corresponding copies of the same molecule type.
%   - A greedy nearest-neighbor match is used within each element block.
%     For the current molecules this should be sufficient.

    p = inputParser;
    addRequired(p, 'X_template', @(x) isnumeric(x) && size(x,2) == 3);
    addRequired(p, 'species_template', @(x) iscell(x) || isstring(x));
    addRequired(p, 'X_target', @(x) isnumeric(x) && size(x,2) == 3);
    addRequired(p, 'species_target', @(x) iscell(x) || isstring(x));
    addParameter(p, 'ReferenceAxis', [0 1 0], @(x) isnumeric(x) && numel(x) == 3);
    addParameter(p, 'DistanceTol', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
    parse(p, X_template, species_template, X_target, species_target, varargin{:});
    opt = p.Results;

    species_template = cellstr(species_template(:));
    species_target = cellstr(species_target(:));

    nT = size(X_template, 1);
    nU = size(X_target, 1);

    if nT ~= nU
        error('builder:match_molecule_atoms_by_frame:AtomCountMismatch', ...
            'Template and target must have the same number of atoms.');
    end

    if numel(species_template) ~= nT || numel(species_target) ~= nU
        error('builder:match_molecule_atoms_by_frame:SpeciesSizeMismatch', ...
            'species_template and species_target must match the number of atoms.');
    end

    if ~same_element_multiset(species_template, species_target)
        error('builder:match_molecule_atoms_by_frame:ElementMismatch', ...
            'Template and target do not contain the same multiset of elements.');
    end

    template_frame = builder.compute_molecule_frame(X_template, ...
        'ReferenceAxis', opt.ReferenceAxis, ...
        'Verbose', false);

    target_frame = builder.compute_molecule_frame(X_target, ...
        'ReferenceAxis', opt.ReferenceAxis, ...
        'Verbose', false);

    Xt = template_frame.local_coords;
    Xu = target_frame.local_coords;

    template_to_target = zeros(nT, 1);
    target_to_template = zeros(nU, 1);
    pair_distance = nan(nT, 1);

    elems = unique(upper(string(species_template)), 'stable');

    for e = 1:numel(elems)
        elem = elems(e);

        idxT = find(strcmpi(species_template, elem));
        idxU = find(strcmpi(species_target, elem));

        if numel(idxT) ~= numel(idxU)
            error('builder:match_molecule_atoms_by_frame:ElementBlockMismatch', ...
                'Element block size mismatch for element %s.', elem);
        end

        DTU = pairwise_distance_matrix(Xt(idxT, :), Xu(idxU, :));

        [mapLocal, distLocal] = greedy_bipartite_match(DTU);

        template_to_target(idxT) = idxU(mapLocal);
        pair_distance(idxT) = distLocal;

        for k = 1:numel(idxT)
            target_to_template(idxU(mapLocal(k))) = idxT(k);
        end
    end

    max_pair_distance = max(pair_distance);
    rms_pair_distance = sqrt(mean(pair_distance.^2));
    is_within_tolerance = all(pair_distance <= opt.DistanceTol);

    match = struct();
    match.template_frame = template_frame;
    match.target_frame = target_frame;
    match.template_to_target = template_to_target;
    match.target_to_template = target_to_template;
    match.template_local_coords = Xt;
    match.target_local_coords = Xu;
    match.pair_distance = pair_distance;
    match.max_pair_distance = max_pair_distance;
    match.rms_pair_distance = rms_pair_distance;
    match.is_within_tolerance = is_within_tolerance;

    if opt.Verbose
        fprintf('Frame-based atom match summary:\n');
        fprintf('  atom count            = %d\n', nT);
        fprintf('  max pair distance     = %.6e\n', max_pair_distance);
        fprintf('  rms pair distance     = %.6e\n', rms_pair_distance);
        fprintf('  within tolerance      = %d\n', is_within_tolerance);
        fprintf('  tolerance             = %.6e\n', opt.DistanceTol);
    end
end


function tf = same_element_multiset(spec1, spec2)
    u1 = unique(upper(string(spec1)));
    u2 = unique(upper(string(spec2)));

    if numel(u1) ~= numel(u2) || ~all(ismember(u1, u2))
        tf = false;
        return;
    end

    tf = true;
    for k = 1:numel(u1)
        tf = tf && (nnz(strcmpi(spec1, u1(k))) == nnz(strcmpi(spec2, u1(k))));
    end
end


function D = pairwise_distance_matrix(X, Y)
    nX = size(X, 1);
    nY = size(Y, 1);
    D = zeros(nX, nY);

    for i = 1:nX
        diff = Y - X(i, :);
        D(i, :) = sqrt(sum(diff.^2, 2)).';
    end
end


function [map, dist] = greedy_bipartite_match(D)
% Map each row to a unique column greedily by smallest remaining distance.
%
% Inputs
%   D     n x n distance matrix
%
% Outputs
%   map   n x 1 column index chosen for each row
%   dist  n x 1 chosen distance for each row

    [nRow, nCol] = size(D);
    if nRow ~= nCol
        error('builder:match_molecule_atoms_by_frame:BadDistanceMatrix', ...
            'Distance matrix must be square for matching.');
    end

    map = zeros(nRow, 1);
    dist = zeros(nRow, 1);

    usedRow = false(nRow, 1);
    usedCol = false(nCol, 1);

    for step = 1:nRow
        bestVal = inf;
        bestI = 0;
        bestJ = 0;

        for i = 1:nRow
            if usedRow(i), continue; end
            for j = 1:nCol
                if usedCol(j), continue; end
                if D(i,j) < bestVal
                    bestVal = D(i,j);
                    bestI = i;
                    bestJ = j;
                end
            end
        end

        if bestI == 0 || bestJ == 0
            error('builder:match_molecule_atoms_by_frame:MatchingFailed', ...
                'Greedy matching failed unexpectedly.');
        end

        map(bestI) = bestJ;
        dist(bestI) = bestVal;
        usedRow(bestI) = true;
        usedCol(bestJ) = true;
    end
end