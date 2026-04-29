function desc = complete_molecule_descriptors_relative_to_reference(sys, refMolID, varargin)
%COMPLETE_MOLECULE_DESCRIPTORS_RELATIVE_TO_REFERENCE
% Build a descriptor table for complete molecules relative to a reference.
%
% desc = builder.complete_molecule_descriptors_relative_to_reference(sys, refMolID)
% desc = builder.complete_molecule_descriptors_relative_to_reference(sys, refMolID, ...)
%
% Optional name-value inputs
%   'StackAxis'        'a' | 'b' | 'c' | numeric 1x3 vector, default 'b'
%   'IncludeReference' logical, default false
%   'IncludeNormals'   logical, default true
%   'Cache'            output of builder.build_complete_molecule_descriptor_cache,
%                      default []
%   'Verbose'          logical, default false
%
% Output desc fields
%   .reference_mol_id
%   .reference_site_indices
%   .reference_com
%   .reference_normal
%   .stack_axis_hat
%   .table
%
% desc.table columns
%   molecule_id
%   is_complete
%   n_sites
%   com_x, com_y, com_z
%   dr_x, dr_y, dr_z
%   d_par
%   d_perp_x, d_perp_y, d_perp_z
%   d_perp
%   distance
%   normal_x, normal_y, normal_z
%   normal_angle_deg

p = inputParser;
addRequired(p, 'sys', @isstruct);
addRequired(p, 'refMolID', @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'StackAxis', 'b', ...
    @(x) ischar(x) || isstring(x) || (isnumeric(x) && numel(x) == 3));
addParameter(p, 'IncludeReference', false, ...
    @(x) islogical(x) && isscalar(x));
addParameter(p, 'IncludeNormals', true, ...
    @(x) islogical(x) && isscalar(x));
addParameter(p, 'Cache', [], ...
    @(x) isempty(x) || isstruct(x));
addParameter(p, 'Verbose', false, ...
    @(x) islogical(x) && isscalar(x));
parse(p, sys, refMolID, varargin{:});

opt = p.Results;

validate_sys(sys);

if isempty(opt.Cache)
    cache = builder.build_complete_molecule_descriptor_cache(sys, ...
        'StackAxis', opt.StackAxis, ...
        'IncludeNormals', opt.IncludeNormals, ...
        'Verbose', false);
else
    cache = opt.Cache;
    validate_cache(cache);

    if opt.IncludeNormals && ~cache.include_normals
        error('builder:complete_molecule_descriptors_relative_to_reference:CacheMissingNormals', ...
            ['IncludeNormals=true was requested, but the supplied descriptor cache ', ...
             'was built with IncludeNormals=false. Rebuild the cache with normals.']);
    end
end

refCacheRow = find(cache.molecule_id(:) == refMolID, 1, 'first');

if isempty(refCacheRow)
    error('builder:complete_molecule_descriptors_relative_to_reference:BadReferenceMolID', ...
        'Reference molecule ID %d is not present among complete molecules.', refMolID);
end

refCOM = cache.com(refCacheRow, :);
stackAxisHat = cache.stack_axis_hat;

if opt.IncludeNormals
    refNormal = cache.normal(refCacheRow, :);
else
    refNormal = [NaN NaN NaN];
end

if opt.IncludeReference
    rows = (1:numel(cache.molecule_id)).';
else
    rows = find(cache.molecule_id(:) ~= refMolID);
end

nMol = numel(rows);

molecule_id = cache.molecule_id(rows);
is_complete = cache.is_complete(rows);
n_sites = cache.n_sites(rows);
com = cache.com(rows, :);

dr = com - refCOM;
d_par = dr * stackAxisHat(:);
d_perp_vec = dr - d_par .* stackAxisHat;
d_perp = vecnorm(d_perp_vec, 2, 2);
distance = vecnorm(dr, 2, 2);

normal = nan(nMol, 3);
normal_angle_deg = nan(nMol, 1);

if opt.IncludeNormals
    normal = cache.normal(rows, :);
    normal_angle_deg = local_unoriented_angle_deg_vectorized(refNormal, normal);
end

[~, order] = sort(distance, 'ascend');

refIdx = builder.site_indices_for_molecule(sys, refMolID);

desc = struct();
desc.reference_mol_id = refMolID;
desc.reference_site_indices = refIdx;
desc.reference_com = refCOM;
desc.reference_normal = refNormal;
desc.stack_axis_hat = stackAxisHat;

desc.table = table( ...
    molecule_id(order), ...
    is_complete(order), ...
    n_sites(order), ...
    com(order,1), com(order,2), com(order,3), ...
    dr(order,1), dr(order,2), dr(order,3), ...
    d_par(order), ...
    d_perp_vec(order,1), d_perp_vec(order,2), d_perp_vec(order,3), ...
    d_perp(order), ...
    distance(order), ...
    normal(order,1), normal(order,2), normal(order,3), ...
    normal_angle_deg(order), ...
    'VariableNames', { ...
        'molecule_id', ...
        'is_complete', ...
        'n_sites', ...
        'com_x', 'com_y', 'com_z', ...
        'dr_x', 'dr_y', 'dr_z', ...
        'd_par', ...
        'd_perp_x', 'd_perp_y', 'd_perp_z', ...
        'd_perp', ...
        'distance', ...
        'normal_x', 'normal_y', 'normal_z', ...
        'normal_angle_deg'});

if opt.Verbose
    fprintf('Descriptor summary:\n');
    fprintf('  reference molecule ID  = %d\n', refMolID);
    fprintf('  complete molecules used= %d\n', nMol);
    fprintf('  stack axis hat         = [%9.4f %9.4f %9.4f]\n', stackAxisHat);

    if opt.IncludeNormals
        fprintf('  reference normal       = [%9.4f %9.4f %9.4f]\n', refNormal);
    end
end

end

% =========================================================================
% Validation / helper functions
% =========================================================================

function validate_sys(sys)

if ~isstruct(sys)
    error('builder:complete_molecule_descriptors_relative_to_reference:BadInput', ...
        'sys must be a struct.');
end

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
    error('builder:complete_molecule_descriptors_relative_to_reference:MissingSitePos', ...
        'sys.site_pos is required and missing/empty.');
end

if ~isfield(sys, 'molecule_table') || isempty(sys.molecule_table)
    error('builder:complete_molecule_descriptors_relative_to_reference:MissingMoleculeTable', ...
        'sys.molecule_table is required and missing/empty.');
end

T = sys.molecule_table;

required = {'molecule_id', 'com', 'n_sites', 'is_complete_in_display'};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(T, name) || isempty(T.(name))
        error('builder:complete_molecule_descriptors_relative_to_reference:BadMoleculeTable', ...
            'sys.molecule_table.%s is required and missing/empty.', name);
    end
end

if size(T.com, 2) ~= 3
    error('builder:complete_molecule_descriptors_relative_to_reference:BadCOM', ...
        'sys.molecule_table.com must be N x 3.');
end

end

function validate_cache(cache)

required = {
    'molecule_id'
    'molecule_table_row'
    'is_complete'
    'n_sites'
    'com'
    'normal'
    'include_normals'
    'stack_axis_hat'
};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(cache, name) || isempty(cache.(name))
        error('builder:complete_molecule_descriptors_relative_to_reference:BadCache', ...
            'cache.%s is required and missing/empty.', name);
    end
end

if size(cache.com, 2) ~= 3
    error('builder:complete_molecule_descriptors_relative_to_reference:BadCacheCOM', ...
        'cache.com must be N x 3.');
end

if size(cache.normal, 2) ~= 3
    error('builder:complete_molecule_descriptors_relative_to_reference:BadCacheNormal', ...
        'cache.normal must be N x 3.');
end

if numel(cache.stack_axis_hat) ~= 3 || norm(cache.stack_axis_hat) == 0
    error('builder:complete_molecule_descriptors_relative_to_reference:BadCacheStackAxis', ...
        'cache.stack_axis_hat must be a nonzero 1x3 vector.');
end

cache.stack_axis_hat = reshape(cache.stack_axis_hat, 1, 3); %#ok<NASGU>

end

function ang = local_unoriented_angle_deg_vectorized(refNormal, normal)

n = size(normal, 1);
ang = nan(n, 1);

if any(isnan(refNormal))
    return;
end

valid = ~any(isnan(normal), 2);

if ~any(valid)
    return;
end

c = abs(normal(valid, :) * refNormal(:));
c = min(max(c, -1), 1);

ang(valid) = acosd(c);

end