function sys = select_active_molecules(sys, activeMolIDs)
%SELECT_ACTIVE_MOLECULES Mark up to two unique molecules as active.
%
% Input
%   sys           working system struct
%   activeMolIDs  vector of unique molecule IDs (length 0, 1, or 2)
%
% Output
%   sys           updated system struct with active-molecule bookkeeping
%
% Added/updated fields
%   sys.active_molecules       unique active molecule IDs
%   sys.site_is_active         N x 1 logical
%   sys.active_site_indices    cell array, one entry per active molecule

if nargin < 2 || isempty(activeMolIDs)
    activeMolIDs = [];
end

activeMolIDs = activeMolIDs(:).';

if numel(activeMolIDs) > 2
    error('At most two active molecules are allowed.');
end

if isempty(sys.site_mol_id)
    error('sys.site_mol_id is missing or empty.');
end

allMolIDs = unique(sys.site_mol_id(:)).';

if ~all(ismember(activeMolIDs, allMolIDs))
    badIDs = activeMolIDs(~ismember(activeMolIDs, allMolIDs));
    error('Active molecule ID(s) not present in system: %s', mat2str(badIDs));
end

sys.active_molecules = activeMolIDs;
sys.site_is_active = false(sys.n_sites, 1);
sys.active_site_indices = cell(numel(activeMolIDs), 1);

for k = 1:numel(activeMolIDs)
    idx = builder.site_indices_for_molecule(sys, activeMolIDs(k));
    sys.site_is_active(idx) = true;
    sys.active_site_indices{k} = idx;
end

end