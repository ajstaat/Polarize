function sys = select_active_molecules(sys, activeMolIDs, varargin)
%SELECT_ACTIVE_MOLECULES Mark selected molecules as active.
%
% sys = builder.select_active_molecules(sys, activeMolIDs)
%
% Active molecules are typically the charged pair. This function does not
% assign charges; it only records active molecule/site bookkeeping.
%
% Output fields set/updated:
%   sys.active_molecules
%   sys.active_site_indices
%   sys.site_is_active

p = inputParser;
addRequired(p, 'sys', @isstruct);
addRequired(p, 'activeMolIDs', @(x) isnumeric(x) && isvector(x));
addParameter(p, 'RequireComplete', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
parse(p, sys, activeMolIDs, varargin{:});

opt = p.Results;

activeMolIDs = unique(activeMolIDs(:), 'stable');

validate_sys(sys);

if opt.RequireComplete
    assert_molecules_complete(sys, activeMolIDs);
end

nSites = size(sys.site_pos, 1);

siteIsActive = false(nSites, 1);
activeSiteIndices = cell(numel(activeMolIDs), 1);

for k = 1:numel(activeMolIDs)
    molID = activeMolIDs(k);
    idx = builder.site_indices_for_molecule(sys, molID);

    if isempty(idx)
        error('builder:select_active_molecules:UnknownMolecule', ...
            'Molecule ID %d was not found.', molID);
    end

    siteIsActive(idx) = true;
    activeSiteIndices{k} = idx(:);

    if opt.Verbose
        fprintf('Active molecule %d: %d sites\n', molID, numel(idx));
    end
end

sys.active_molecules = activeMolIDs(:);
sys.active_site_indices = activeSiteIndices;
sys.site_is_active = siteIsActive;

end

function validate_sys(sys)

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
    error('builder:select_active_molecules:MissingSitePos', ...
        'sys.site_pos is required.');
end

n = size(sys.site_pos, 1);

if ~isfield(sys, 'site_mol_id') || numel(sys.site_mol_id) ~= n
    error('builder:select_active_molecules:BadSiteMolID', ...
        'sys.site_mol_id must have one entry per site.');
end

if ~isfield(sys, 'molecule_table') || isempty(sys.molecule_table) || ...
        ~isfield(sys.molecule_table, 'molecule_id')
    error('builder:select_active_molecules:MissingMoleculeTable', ...
        'sys.molecule_table.molecule_id is required.');
end

end

function assert_molecules_complete(sys, molIDs)

T = sys.molecule_table;

if ~isfield(T, 'is_complete_in_display')
    error('builder:select_active_molecules:MissingCompletenessFlag', ...
        'sys.molecule_table.is_complete_in_display is required.');
end

for k = 1:numel(molIDs)
    row = find(T.molecule_id == molIDs(k), 1, 'first');

    if isempty(row)
        error('builder:select_active_molecules:UnknownMolecule', ...
            'Molecule ID %d was not found.', molIDs(k));
    end

    if ~T.is_complete_in_display(row)
        error('builder:select_active_molecules:IncompleteMolecule', ...
            'Molecule ID %d is not complete in the displayed supercell.', molIDs(k));
    end
end

end