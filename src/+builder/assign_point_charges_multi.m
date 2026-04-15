function sys = assign_point_charges_multi(sys, chargePatterns)
%ASSIGN_POINT_CHARGES_MULTI Assign distinct charge patterns to active molecules.
%
% Inputs
%   sys            working system struct
%   chargePatterns cell array of structs, each with fields:
%                    .mol_id     unique molecule ID in the built supercell
%                    .site_label cell array of site labels within that molecule
%                    .delta_q    numeric vector of charge changes
%
% Behavior
%   - Charges are assigned only to the molecule named by each pattern.mol_id
%   - Matching is done by site label within that molecule
%   - Sites not mentioned in any pattern get zero charge
%   - A given molecule may appear at most once in chargePatterns
%
% Added/updated fields
%   sys.site_charge
%   sys.charge_patterns
%
% Example
%   chargePatterns = {
%       struct('mol_id', 6, 'site_label', {{'C1','N1'}}, 'delta_q', [ 0.10, -0.10]), ...
%       struct('mol_id', 11,'site_label', {{'C1','N1'}}, 'delta_q', [-0.10,  0.10])
%   };

    if nargin < 2 || isempty(chargePatterns)
        error('chargePatterns must be provided and non-empty.');
    end

    if ~iscell(chargePatterns)
        error('chargePatterns must be a cell array of structs.');
    end

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

    activeIDs = sys.active_molecules(:).';
    targetMolIDs = zeros(numel(chargePatterns), 1);

    for p = 1:numel(chargePatterns)
        pat = chargePatterns{p};

        if ~isstruct(pat)
            error('Each entry of chargePatterns must be a struct.');
        end
        if ~isfield(pat, 'mol_id') || isempty(pat.mol_id)
            error('Each charge pattern must include a scalar mol_id.');
        end
        if ~isscalar(pat.mol_id) || ~isnumeric(pat.mol_id)
            error('charge pattern mol_id must be a numeric scalar.');
        end

        targetMolIDs(p) = pat.mol_id;

        if ~ismember(pat.mol_id, activeIDs)
            error('Charge pattern mol_id=%d is not in sys.active_molecules.', pat.mol_id);
        end

        builder.validate_charge_pattern(pat);

        idxMol = builder.site_indices_for_molecule(sys, pat.mol_id);
        molLabels = sys.site_label(idxMol);

        patternLabels = pat.site_label(:);
        patternDQ = pat.delta_q(:);

        for j = 1:numel(patternLabels)
            thisLabel = patternLabels{j};
            thisDQ = patternDQ(j);

            localMatch = strcmp(molLabels, thisLabel);
            nMatch = sum(localMatch);

            if nMatch == 0
                error(['Charge pattern label "%s" not found in molecule %d. ', ...
                       'Check site labels or molecule definition.'], ...
                       thisLabel, pat.mol_id);
            elseif nMatch > 1
                error(['Charge pattern label "%s" appears multiple times in molecule %d. ', ...
                       'Site labels must uniquely identify sites within a molecule.'], ...
                       thisLabel, pat.mol_id);
            end

            globalIdx = idxMol(localMatch);
            sys.site_charge(globalIdx) = thisDQ;
        end
    end

    if numel(unique(targetMolIDs)) ~= numel(targetMolIDs)
        error('Each molecule may appear at most once in chargePatterns.');
    end

    sys.charge_patterns = chargePatterns;
end