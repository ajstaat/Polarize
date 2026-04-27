function sys = remove_molecules(sys, removeMolIDs, varargin)
%REMOVE_MOLECULES Remove molecule components from a system.
%
% sys = builder.remove_molecules(sys, removeMolIDs)
%
% This removes all sites whose sys.site_mol_id is in removeMolIDs and
% updates molecule_table.site_indices to the new site indexing.
%
% Molecule IDs are preserved for remaining molecules; they are not compacted.
% This is intentional because downstream code may refer to selected molecule
% IDs from the pre-removal system.

p = inputParser;
addRequired(p, 'sys', @isstruct);
addRequired(p, 'removeMolIDs', @(x) isnumeric(x) && isvector(x));
addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
parse(p, sys, removeMolIDs, varargin{:});

removeMolIDs = unique(removeMolIDs(:), 'stable');

validate_sys(sys);

nOld = size(sys.site_pos, 1);

removeMask = ismember(sys.site_mol_id(:), removeMolIDs);
keepMask = ~removeMask;

if ~any(removeMask)
    if p.Results.Verbose
        fprintf('remove_molecules: no matching molecules to remove.\n');
    end
    return;
end

oldToNew = zeros(nOld, 1);
oldToNew(keepMask) = 1:nnz(keepMask);

siteFields = local_site_fields(sys, nOld);

for k = 1:numel(siteFields)
    name = siteFields{k};
    sys.(name) = local_slice_site_field(sys.(name), keepMask);
end

% Square site-site matrices, if present.
matrixFields = {
    'supercell_pbc_bond_graph'
};

for k = 1:numel(matrixFields)
    name = matrixFields{k};

    if isfield(sys, name) && ~isempty(sys.(name)) && ...
            isequal(size(sys.(name)), [nOld nOld])
        sys.(name) = sys.(name)(keepMask, keepMask);
    end
end

sys.n_sites = nnz(keepMask);

% Update molecule table.
if isfield(sys, 'molecule_table') && ~isempty(sys.molecule_table)
    sys.molecule_table = local_update_molecule_table(sys.molecule_table, removeMolIDs, oldToNew);
end

% Update active bookkeeping if present.
if isfield(sys, 'site_is_active') && numel(sys.site_is_active) == sys.n_sites
    activeIDs = unique(sys.site_mol_id(logical(sys.site_is_active)), 'stable');
    sys.active_molecules = activeIDs(:);

    sys.active_site_indices = cell(numel(activeIDs), 1);
    for k = 1:numel(activeIDs)
        sys.active_site_indices{k} = find(sys.site_mol_id == activeIDs(k));
    end
else
    sys.site_is_active = false(sys.n_sites, 1);
    sys.active_molecules = [];
    sys.active_site_indices = {};
end

if isfield(sys, 'removed_molecules') && ~isempty(sys.removed_molecules)
    sys.removed_molecules = unique([sys.removed_molecules(:); removeMolIDs(:)], 'stable');
else
    sys.removed_molecules = removeMolIDs(:);
end

if p.Results.Verbose
    fprintf('remove_molecules: removed %d molecule(s), %d site(s).\n', ...
        numel(removeMolIDs), nnz(removeMask));
    fprintf('remove_molecules: remaining sites = %d\n', sys.n_sites);
end

end

% =========================================================================
% Helpers
% =========================================================================

function validate_sys(sys)

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
    error('builder:remove_molecules:MissingSitePos', ...
        'sys.site_pos is required.');
end

n = size(sys.site_pos, 1);

if ~isfield(sys, 'site_mol_id') || numel(sys.site_mol_id) ~= n
    error('builder:remove_molecules:BadSiteMolID', ...
        'sys.site_mol_id must have one entry per site.');
end

end

function fields = local_site_fields(sys, nSites)

allNames = fieldnames(sys);
fields = {};

for k = 1:numel(allNames)
    name = allNames{k};
    value = sys.(name);

    if isstruct(value) || istable(value)
        continue;
    end

    if isvector(value) && numel(value) == nSites
        fields{end+1} = name; %#ok<AGROW>
    elseif isnumeric(value) || islogical(value)
        if size(value, 1) == nSites && ~isequal(size(value), [nSites nSites])
            fields{end+1} = name; %#ok<AGROW>
        end
    elseif iscell(value)
        if size(value, 1) == nSites || numel(value) == nSites
            fields{end+1} = name; %#ok<AGROW>
        end
    end
end

fields = unique(fields, 'stable');

end

function y = local_slice_site_field(x, keepMask)

if isvector(x) && numel(x) == numel(keepMask)
    y = x(keepMask);
elseif size(x, 1) == numel(keepMask)
    y = x(keepMask, :);
else
    y = x;
end

end

function T = local_update_molecule_table(T, removeMolIDs, oldToNew)

if ~isfield(T, 'molecule_id')
    return;
end

keepRows = ~ismember(T.molecule_id(:), removeMolIDs(:));

names = fieldnames(T);
nMolOld = numel(T.molecule_id);

for k = 1:numel(names)
    name = names{k};
    value = T.(name);

    if isvector(value) && numel(value) == nMolOld
        T.(name) = value(keepRows);
    elseif size(value, 1) == nMolOld
        T.(name) = value(keepRows, :);
    end
end

if isfield(T, 'site_indices')
    for r = 1:numel(T.site_indices)
        oldIdx = T.site_indices{r};
        newIdx = oldToNew(oldIdx);
        newIdx = newIdx(newIdx > 0);
        T.site_indices{r} = newIdx(:);
    end
end

end