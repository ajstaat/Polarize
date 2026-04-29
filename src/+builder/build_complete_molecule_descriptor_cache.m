function cache = build_complete_molecule_descriptor_cache(sys, varargin)
%BUILD_COMPLETE_MOLECULE_DESCRIPTOR_CACHE Cache complete-molecule descriptors.
%
% cache = builder.build_complete_molecule_descriptor_cache(sys, ...)
%
% Builds molecule-level data used by neighbor selectors:
%   molecule IDs
%   molecule_table rows
%   COMs
%   n_sites
%   optional best-fit molecular normals
%   stack-axis unit vector
%
% The point is to compute these quantities once per system, then reuse them
% for each reference molecule in:
%
%   builder.complete_molecule_descriptors_relative_to_reference
%
% Options
%   'StackAxis'       'a' | 'b' | 'c' | numeric 1x3 vector, default 'b'
%   'IncludeNormals'  logical, default true
%   'Verbose'         logical, default false

p = inputParser;
addRequired(p, 'sys', @isstruct);
addParameter(p, 'StackAxis', 'b', ...
    @(x) ischar(x) || isstring(x) || (isnumeric(x) && numel(x) == 3));
addParameter(p, 'IncludeNormals', true, ...
    @(x) islogical(x) && isscalar(x));
addParameter(p, 'Verbose', false, ...
    @(x) islogical(x) && isscalar(x));
parse(p, sys, varargin{:});

opt = p.Results;

validate_sys(sys);

T = sys.molecule_table;

completeRows = find(logical(T.is_complete_in_display(:)));

if isempty(completeRows)
    error('builder:build_complete_molecule_descriptor_cache:NoCompleteMolecules', ...
        'No complete molecules are available in sys.molecule_table.');
end

moleculeID = T.molecule_id(completeRows);
com = T.com(completeRows, :);
nSites = T.n_sites(completeRows);
isComplete = logical(T.is_complete_in_display(completeRows));

stackAxisHat = local_resolve_stack_axis(sys, opt.StackAxis);

nMol = numel(moleculeID);

normal = nan(nMol, 3);

if opt.IncludeNormals
    for k = 1:nMol
        idx = builder.site_indices_for_molecule(sys, moleculeID(k));
        normal(k, :) = local_estimate_molecule_normal(sys.site_pos(idx, :));
    end
end

cache = struct();
cache.molecule_id = moleculeID(:);
cache.molecule_table_row = completeRows(:);
cache.is_complete = isComplete(:);
cache.n_sites = nSites(:);
cache.com = com;
cache.normal = normal;
cache.include_normals = opt.IncludeNormals;
cache.stack_axis_hat = stackAxisHat;
cache.stack_axis_input = opt.StackAxis;

if opt.Verbose
    fprintf('Complete-molecule descriptor cache:\n');
    fprintf('  complete molecules = %d\n', nMol);
    fprintf('  include normals    = %d\n', opt.IncludeNormals);
    fprintf('  stack axis hat     = [%9.4f %9.4f %9.4f]\n', stackAxisHat);
end

end

% =========================================================================
% Validation / geometry helpers
% =========================================================================

function validate_sys(sys)

if ~isstruct(sys)
    error('builder:build_complete_molecule_descriptor_cache:BadInput', ...
        'sys must be a struct.');
end

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
    error('builder:build_complete_molecule_descriptor_cache:MissingSitePos', ...
        'sys.site_pos is required and missing/empty.');
end

if ~isfield(sys, 'molecule_table') || isempty(sys.molecule_table)
    error('builder:build_complete_molecule_descriptor_cache:MissingMoleculeTable', ...
        'sys.molecule_table is required and missing/empty.');
end

T = sys.molecule_table;

required = {'molecule_id', 'com', 'n_sites', 'is_complete_in_display'};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(T, name) || isempty(T.(name))
        error('builder:build_complete_molecule_descriptor_cache:BadMoleculeTable', ...
            'sys.molecule_table.%s is required and missing/empty.', name);
    end
end

if size(T.com, 2) ~= 3
    error('builder:build_complete_molecule_descriptor_cache:BadCOM', ...
        'sys.molecule_table.com must be N x 3.');
end

end

function stackAxisHat = local_resolve_stack_axis(sys, axisSpec)

if ischar(axisSpec) || (isstring(axisSpec) && isscalar(axisSpec))
    axisSpec = lower(char(string(axisSpec)));

    if ~isfield(sys, 'super_lattice') || isempty(sys.super_lattice)
        error('builder:build_complete_molecule_descriptor_cache:MissingSuperLattice', ...
            'sys.super_lattice is required for symbolic stack-axis selection.');
    end

    switch axisSpec
        case 'a'
            v = sys.super_lattice(1, :);
        case 'b'
            v = sys.super_lattice(2, :);
        case 'c'
            v = sys.super_lattice(3, :);
        otherwise
            error('builder:build_complete_molecule_descriptor_cache:BadStackAxis', ...
                'StackAxis must be ''a'', ''b'', ''c'', or a numeric 1x3 vector.');
    end
elseif isnumeric(axisSpec) && numel(axisSpec) == 3
    v = reshape(axisSpec, 1, 3);
else
    error('builder:build_complete_molecule_descriptor_cache:BadStackAxis', ...
        'StackAxis must be ''a'', ''b'', ''c'', or a numeric 1x3 vector.');
end

nv = norm(v);

if nv == 0
    error('builder:build_complete_molecule_descriptor_cache:ZeroStackAxis', ...
        'Resolved stack axis has zero norm.');
end

stackAxisHat = v / nv;

end

function nHat = local_estimate_molecule_normal(X)
% Best-fit plane normal from site positions X (N x 3).

if size(X, 1) < 3
    nHat = [NaN NaN NaN];
    return;
end

Xc = X - mean(X, 1);

[~, S, V] = svd(Xc, 'econ');

sing = diag(S);

if numel(sing) < 3
    nHat = [NaN NaN NaN];
    return;
end

nHat = V(:, end).';

nn = norm(nHat);

if nn == 0
    nHat = [NaN NaN NaN];
else
    nHat = nHat / nn;
end

end