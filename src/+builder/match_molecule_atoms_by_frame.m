function map = match_molecule_atoms_by_frame(template, target, varargin)
%MATCH_MOLECULE_ATOMS_BY_FRAME Match template atoms to target atoms by oriented frame.
%
% map = builder.match_molecule_atoms_by_frame(template, target)
% map = builder.match_molecule_atoms_by_frame(template, target, Name, Value)
%
% Inputs
%   template struct with fields:
%       .site_pos   N x 3 Cartesian coordinates
%       .site_type  N x 1 atom/species labels
%
%   target struct with fields:
%       .site_pos   N x 3 Cartesian coordinates
%       .site_type  N x 1 atom/species labels
%
% Options
%   'DistanceTol'
%       maximum allowed same-element local-coordinate match distance
%       default: 0.2
%
%   'ReferenceAxis'
%       crystal/reference direction used to fix molecular-normal sign e3
%       default: [0 0 1]
%
%   'PrimaryAxis'
%       crystal/reference direction projected into the molecular plane to
%       define positive e1
%       default: [1 0 0]
%
%   'AmbiguityTol'
%       if two assignments of same-element atoms have objective values
%       within this tolerance, error rather than silently choosing
%       default: 1e-10
%
% Output map fields
%   .template_to_target
%       N x 1 vector; target index for each template atom
%
%   .target_to_template
%       N x 1 vector; template index for each target atom
%
%   .distance
%       N x 1 matched local-coordinate distances
%
%   .max_distance
%   .rms_distance
%   .template_frame
%   .target_frame
%   .template_local
%   .target_local
%
% Notes
%   Coordinates must already be in the same length unit. In the current
%   builder workflow, sys.site_pos is bohr, so template coordinates used for
%   charge mapping should also be bohr.
%
%   This matcher preserves orientation. It is appropriate for charge-file
%   mapping where symmetry-equivalent atoms may still need to be
%   distinguished by molecular face / crystal-axis orientation.

p = inputParser;
addRequired(p, 'template', @isstruct);
addRequired(p, 'target', @isstruct);
addParameter(p, 'DistanceTol', 0.2, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'ReferenceAxis', [0 0 1], ...
    @(x) isnumeric(x) && numel(x) == 3 && norm(x) > 0);
addParameter(p, 'PrimaryAxis', [1 0 0], ...
    @(x) isnumeric(x) && numel(x) == 3 && norm(x) > 0);
addParameter(p, 'AmbiguityTol', 1e-10, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
parse(p, template, target, varargin{:});

opt = p.Results;

template = local_normalize_molecule(template, 'template');
target = local_normalize_molecule(target, 'target');

n = size(template.site_pos, 1);

if size(target.site_pos, 1) ~= n
    error('builder:match_molecule_atoms_by_frame:SizeMismatch', ...
        'Template and target must contain the same number of atoms.');
end

local_assert_same_composition(template.site_type, target.site_type);

refAxis = reshape(opt.ReferenceAxis, 1, 3);
refAxis = refAxis / norm(refAxis);

primaryAxis = reshape(opt.PrimaryAxis, 1, 3);
primaryAxis = primaryAxis / norm(primaryAxis);

templateFrame = local_compute_oriented_frame(template.site_pos, refAxis, primaryAxis);
targetFrame = local_compute_oriented_frame(target.site_pos, refAxis, primaryAxis);

templateLocal = local_to_frame_coordinates(template.site_pos, templateFrame);
targetLocal = local_to_frame_coordinates(target.site_pos, targetFrame);

candidate = local_match_by_type_and_local_coordinates( ...
    templateLocal, template.site_type, ...
    targetLocal, target.site_type, ...
    opt.DistanceTol, opt.AmbiguityTol);

map = struct();

map.template_to_target = candidate.template_to_target(:);
map.target_to_template = candidate.target_to_template(:);

map.distance = candidate.distance(:);
map.max_distance = candidate.max_distance;
map.rms_distance = candidate.rms_distance;

map.template_frame = templateFrame;
map.target_frame = targetFrame;

map.template_local = templateLocal;
map.target_local = targetLocal;

end

% =========================================================================
% Molecule/frame helpers
% =========================================================================

function mol = local_normalize_molecule(mol, role)

required = {'site_pos', 'site_type'};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(mol, name) || isempty(mol.(name))
        error('builder:match_molecule_atoms_by_frame:MissingField', ...
            '%s.%s is required and missing/empty.', role, name);
    end
end

if ~isnumeric(mol.site_pos) || size(mol.site_pos, 2) ~= 3
    error('builder:match_molecule_atoms_by_frame:BadSitePos', ...
        '%s.site_pos must be N x 3 numeric.', role);
end

mol.site_type = local_to_cell_column(mol.site_type);

if numel(mol.site_type) ~= size(mol.site_pos, 1)
    error('builder:match_molecule_atoms_by_frame:BadSiteType', ...
        '%s.site_type must have one entry per atom.', role);
end

end

function frame = local_compute_oriented_frame(X, referenceAxis, primaryAxis)
%LOCAL_COMPUTE_ORIENTED_FRAME Deterministic orientable molecular frame.
%
% e3:
%   best-fit molecular normal, signed by referenceAxis
%
% e1:
%   primaryAxis projected into molecular plane
%
% e2:
%   cross(e3, e1)

com = mean(X, 1);
X0 = X - com;

[~, ~, V] = svd(X0, 0);

% For approximately planar molecules, smallest-variance direction is normal.
e3 = V(:, end).';

if dot(e3, referenceAxis) < 0
    e3 = -e3;
end

e3 = e3 / norm(e3);

% Define e1 from projected primary axis. This is the key orientability
% convention relative to the cell/crystal.
e1 = primaryAxis - dot(primaryAxis, e3) * e3;

if norm(e1) < 1e-10
    % Fallback: use largest-variance principal direction if primaryAxis is
    % too close to the molecular normal.
    e1 = V(:, 1).';
    e1 = e1 - dot(e1, e3) * e3;

    if dot(e1, primaryAxis) < 0
        e1 = -e1;
    end
end

e1 = e1 / norm(e1);

e2 = cross(e3, e1);
e2 = e2 / norm(e2);

% Re-orthogonalize e1 to avoid tiny numerical drift.
e1 = cross(e2, e3);
e1 = e1 / norm(e1);

frame = struct();
frame.com = com;
frame.e1 = e1;
frame.e2 = e2;
frame.e3 = e3;

R = [
    e1
    e2
    e3
];

frame.local_coords = X0 * R.';

end

function Xlocal = local_to_frame_coordinates(X, frame)

R = [
    frame.e1(:).'
    frame.e2(:).'
    frame.e3(:).'
];

X0 = X - frame.com;
Xlocal = X0 * R.';

end

% =========================================================================
% Assignment helpers
% =========================================================================

function candidate = local_match_by_type_and_local_coordinates(templateLocal, templateTypes, ...
    targetLocal, targetTypes, distanceTol, ambiguityTol)

n = size(templateLocal, 1);

templateToTarget = zeros(n, 1);
targetToTemplate = zeros(n, 1);
distance = zeros(n, 1);

types = unique(templateTypes, 'stable');

for t = 1:numel(types)
    typ = types{t};

    iT = find(strcmp(templateTypes, typ));
    iU = find(strcmp(targetTypes, typ));

    if numel(iT) ~= numel(iU)
        error('builder:match_molecule_atoms_by_frame:CompositionMismatch', ...
            'Type counts do not match for "%s".', typ);
    end

    D = zeros(numel(iT), numel(iU));

    for i = 1:numel(iT)
        for j = 1:numel(iU)
            D(i, j) = norm(templateLocal(iT(i), :) - targetLocal(iU(j), :));
        end
    end

    assignment = local_best_assignment(D, ambiguityTol);

    for row = 1:size(assignment.pairs, 1)
        iLocal = assignment.pairs(row, 1);
        jLocal = assignment.pairs(row, 2);

        i = iT(iLocal);
        j = iU(jLocal);

        d = D(iLocal, jLocal);

        if d > distanceTol
            error('builder:match_molecule_atoms_by_frame:DistanceTooLarge', ...
                ['Best match for template atom %d type "%s" has local-frame distance %.6g, ' ...
                 'which exceeds DistanceTol %.6g.'], ...
                i, typ, d, distanceTol);
        end

        templateToTarget(i) = j;
        targetToTemplate(j) = i;
        distance(i) = d;
    end
end

if any(templateToTarget == 0) || any(targetToTemplate == 0)
    error('builder:match_molecule_atoms_by_frame:IncompleteAssignment', ...
        'Atom assignment did not cover all atoms.');
end

candidate = struct();
candidate.template_to_target = templateToTarget;
candidate.target_to_template = targetToTemplate;
candidate.distance = distance;
candidate.max_distance = max(distance);
candidate.rms_distance = sqrt(mean(distance.^2));

end

function assignment = local_best_assignment(D, ambiguityTol)
%LOCAL_BEST_ASSIGNMENT Find minimum sum-of-squares assignment.
%
% For small same-element groups, use exhaustive permutations. This avoids
% greedy failures in symmetric molecules. For larger groups, fall back to
% greedy with a warning.

[nRows, nCols] = size(D);

if nRows ~= nCols
    error('builder:match_molecule_atoms_by_frame:AssignmentShape', ...
        'Assignment matrix must be square.');
end

if nRows == 1
    assignment = struct();
    assignment.pairs = [1 1];
    assignment.score = D(1,1)^2;
    assignment.second_score = inf;
    return;
end

if nRows <= 8
    P = perms(1:nRows);
    nPerm = size(P, 1);

    scores = zeros(nPerm, 1);

    for p = 1:nPerm
        s = 0;

        for i = 1:nRows
            s = s + D(i, P(p, i))^2;
        end

        scores(p) = s;
    end

    [scoresSorted, order] = sort(scores, 'ascend');

    bestIdx = order(1);
    bestPerm = P(bestIdx, :);

    if numel(scoresSorted) > 1
        secondScore = scoresSorted(2);
    else
        secondScore = inf;
    end

    if isfinite(secondScore) && abs(secondScore - scoresSorted(1)) <= ambiguityTol
        error('builder:match_molecule_atoms_by_frame:AmbiguousAssignment', ...
            ['Frame-based atom mapping is ambiguous. Best assignment score %.6g, ' ...
             'second-best %.6g.'], scoresSorted(1), secondScore);
    end

    pairs = zeros(nRows, 2);
    for i = 1:nRows
        pairs(i, :) = [i bestPerm(i)];
    end

    assignment = struct();
    assignment.pairs = pairs;
    assignment.score = scoresSorted(1);
    assignment.second_score = secondScore;

else
    warning('builder:match_molecule_atoms_by_frame:GreedyFallback', ...
        'Same-element group has %d atoms; using greedy assignment fallback.', nRows);

    pairs = local_greedy_assignment(D);

    assignment = struct();
    assignment.pairs = pairs;
    assignment.score = NaN;
    assignment.second_score = NaN;
end

end

function pairs = local_greedy_assignment(D)

[nRows, nCols] = size(D);

pairs = zeros(nRows, 2);

usedRows = false(nRows, 1);
usedCols = false(nCols, 1);

for k = 1:nRows
    Dwork = D;
    Dwork(usedRows, :) = inf;
    Dwork(:, usedCols) = inf;

    [~, linearIdx] = min(Dwork(:));
    [i, j] = ind2sub(size(Dwork), linearIdx);

    if ~isfinite(Dwork(i, j))
        error('builder:match_molecule_atoms_by_frame:AssignmentFailed', ...
            'Could not complete assignment.');
    end

    pairs(k, :) = [i j];

    usedRows(i) = true;
    usedCols(j) = true;
end

end

% =========================================================================
% Text/composition helpers
% =========================================================================

function c = local_to_cell_column(x)

if isstring(x)
    c = cellstr(x(:));
elseif iscellstr(x)
    c = x(:);
elseif iscell(x)
    c = x(:);
else
    error('builder:match_molecule_atoms_by_frame:BadTextMetadata', ...
        'site_type must be string/cellstr/cell.');
end

c = cellfun(@(s) char(string(s)), c, 'UniformOutput', false);
c = c(:);

end

function local_assert_same_composition(typesA, typesB)

uA = unique(typesA, 'stable');
uB = unique(typesB, 'stable');

if numel(uA) ~= numel(uB) || ~all(ismember(uA, uB)) || ~all(ismember(uB, uA))
    error('builder:match_molecule_atoms_by_frame:CompositionMismatch', ...
        'Template and target have different element/type sets.');
end

for k = 1:numel(uA)
    typ = uA{k};

    if nnz(strcmp(typesA, typ)) ~= nnz(strcmp(typesB, typ))
        error('builder:match_molecule_atoms_by_frame:CompositionMismatch', ...
            'Template and target have different counts for type "%s".', typ);
    end
end

end