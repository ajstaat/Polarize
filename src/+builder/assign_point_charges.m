function sys = assign_point_charges(sys, chargePattern)
%ASSIGN_POINT_CHARGES Assign a charge-change pattern to active molecules.
%
% Inputs
%   sys            working system struct
%   chargePattern  struct with fields:
%                  .site_label   cell array of site labels in one molecule
%                  .delta_q      matching numeric vector of charge changes
%
% Behavior
%   - Charges are assigned only to sites in sys.active_molecules
%   - Matching is done by site label within each active molecule
%   - Sites not mentioned in the pattern get zero charge
%
% Added/updated fields
%   sys.site_charge
%   sys.charge_pattern
%
% Example
%   chargePattern.site_label = {'C1','N1'};
%   chargePattern.delta_q = [0.10, -0.10];

builder.validate_charge_pattern(chargePattern);

if ~isfield(sys, 'active_molecules') || isempty(sys.active_molecules)
    error('No active molecules are selected. Call builder.select_active_molecules first.');
end

if ~isfield(sys, 'site_label') || isempty(sys.site_label)
    error('sys.site_label is required to assign charges by site label.');
end

if ~isfield(sys, 'site_mol_id') || isempty(sys.site_mol_id)
    error('sys.site_mol_id is missing or empty.');
end

sys = builder.clear_site_charges(sys);

patternLabels = chargePattern.site_label(:);
patternDQ = chargePattern.delta_q(:);

for k = 1:numel(sys.active_molecules)
    molID = sys.active_molecules(k);
    idxMol = builder.site_indices_for_molecule(sys, molID);

    molLabels = sys.site_label(idxMol);

    for j = 1:numel(patternLabels)
        thisLabel = patternLabels{j};
        thisDQ = patternDQ(j);

        localMatch = strcmp(molLabels, thisLabel);
        nMatch = sum(localMatch);

        if nMatch == 0
            error(['Charge pattern label "%s" not found in active molecule %d. ', ...
                   'Check site labels or molecule definition.'], thisLabel, molID);
        elseif nMatch > 1
            error(['Charge pattern label "%s" appears multiple times in active molecule %d. ', ...
                   'Site labels must uniquely identify sites within a molecule.'], thisLabel, molID);
        end

        globalIdx = idxMol(localMatch);
        sys.site_charge(globalIdx) = thisDQ;
    end
end

sys.charge_pattern = chargePattern;

end