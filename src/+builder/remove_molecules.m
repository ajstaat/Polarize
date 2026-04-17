function sys = remove_molecules(sys, removeMolIDs)
%REMOVE_MOLECULES Remove selected unique molecules from the working system.
%
% Input
%   sys           working system struct
%   removeMolIDs  vector of unique molecule IDs to remove
%
% Output
%   sys           updated system struct with selected molecules removed
%
% Notes
%   - Removal is done on the current working arrays.
%   - This function assumes sys.site_mol_id contains unique molecule IDs.
%   - Active molecules cannot also be removed.

    if nargin < 2 || isempty(removeMolIDs)
        removeMolIDs = [];
    end

    removeMolIDs = unique(removeMolIDs(:).', 'stable');

    if isempty(removeMolIDs)
        return;
    end

    if ~isfield(sys, 'site_mol_id') || isempty(sys.site_mol_id)
        error('builder:remove_molecules:MissingSiteMolId', ...
            'sys.site_mol_id is missing or empty.');
    end

    allMolIDs = unique(sys.site_mol_id(:)).';

    if ~all(ismember(removeMolIDs, allMolIDs))
        badIDs = removeMolIDs(~ismember(removeMolIDs, allMolIDs));
        error('builder:remove_molecules:BadMolId', ...
            'Remove molecule ID(s) not present in system: %s', mat2str(badIDs));
    end

    if isfield(sys, 'active_molecules') && ~isempty(sys.active_molecules)
        overlap = intersect(removeMolIDs, sys.active_molecules);
        if ~isempty(overlap)
            error('builder:remove_molecules:ActiveOverlap', ...
                'Cannot remove molecule(s) that are also active: %s', mat2str(overlap));
        end
    end

    keepMask = ~ismember(sys.site_mol_id, removeMolIDs);

    % Apply mask to site-level numeric/logical arrays
    siteFields = { ...
        'site_pos', ...
        'site_frac', ...
        'site_mol_id', ...
        'base_mol_id', ...
        'unique_mol_id', ...
        'site_is_active', ...
        'site_is_polarizable', ...
        'site_alpha', ...
        'site_charge', ...
        'cell_shift', ...
        'image_id', ...
        'unit_site_index'};

    for f = 1:numel(siteFields)
        name = siteFields{f};
        if isfield(sys, name) && ~isempty(sys.(name))
            value = sys.(name);

            if isvector(value) && numel(value) == numel(keepMask)
                sys.(name) = value(keepMask, :);
            elseif size(value, 1) == numel(keepMask)
                sys.(name) = value(keepMask, :);
            end
        end
    end

    % Apply mask to site-level cell arrays
    cellFields = {'site_label', 'site_type'};
    for f = 1:numel(cellFields)
        name = cellFields{f};
        if isfield(sys, name) && ~isempty(sys.(name))
            value = sys.(name);
            if numel(value) == numel(keepMask)
                sys.(name) = value(keepMask);
            end
        end
    end

    % Track removed molecules
    if ~isfield(sys, 'removed_molecules') || isempty(sys.removed_molecules)
        sys.removed_molecules = removeMolIDs;
    else
        sys.removed_molecules = unique([sys.removed_molecules(:).', removeMolIDs], 'stable');
    end

    % Update site count
    sys.n_sites = sum(keepMask);

    % Rebuild active-site index cache if needed
    if isfield(sys, 'active_molecules') && ~isempty(sys.active_molecules)
        sys.active_site_indices = cell(numel(sys.active_molecules), 1);
        for k = 1:numel(sys.active_molecules)
            sys.active_site_indices{k} = builder.site_indices_for_molecule(sys, sys.active_molecules(k));
        end
    else
        sys.active_site_indices = {};
        if isfield(sys, 'site_is_active')
            sys.site_is_active = false(sys.n_sites, 1);
        end
    end

    % Trim molecule_table without renumbering unique molecule IDs
    if isfield(sys, 'molecule_table') && isstruct(sys.molecule_table) && ...
       isfield(sys.molecule_table, 'unique_mol_id')

        keepMol = ~ismember(sys.molecule_table.unique_mol_id, removeMolIDs);

        sys.molecule_table.unique_mol_id = sys.molecule_table.unique_mol_id(keepMol);
        sys.molecule_table.base_mol_id   = sys.molecule_table.base_mol_id(keepMol);
        sys.molecule_table.cell_shift    = sys.molecule_table.cell_shift(keepMol, :);
        sys.molecule_table.site_indices  = sys.molecule_table.site_indices(keepMol);

        % Recompute site_indices because site numbering changed after masking
        for m = 1:numel(sys.molecule_table.unique_mol_id)
            thisMol = sys.molecule_table.unique_mol_id(m);
            sys.molecule_table.site_indices{m} = find(sys.site_mol_id == thisMol);
        end
    end
end