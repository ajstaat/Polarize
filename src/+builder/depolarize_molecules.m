function sys = depolarize_molecules(sys, molIDs)
%DEPOLARIZE_MOLECULES Zero site polarizability on selected unique molecules.
%
% Inputs
%   sys    working system struct
%   molIDs vector of unique molecule IDs
%
% Effects
%   - sets sys.site_alpha(idx) = 0
%   - sets sys.site_is_polarizable(idx) = false
%
% Notes
%   molIDs are unique molecule IDs in the built supercell, i.e. values
%   compatible with sys.site_mol_id and builder.select_active_molecules.

    if nargin < 2 || isempty(molIDs)
        return;
    end

    if ~isfield(sys, 'site_mol_id') || isempty(sys.site_mol_id)
        error('sys.site_mol_id is missing.');
    end
    if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
        error('sys.site_alpha is missing.');
    end
    if ~isfield(sys, 'site_is_polarizable') || isempty(sys.site_is_polarizable)
        error('sys.site_is_polarizable is missing.');
    end

    molIDs = molIDs(:);
    idx = ismember(sys.site_mol_id, molIDs);

    sys.site_alpha(idx) = 0;
    sys.site_is_polarizable(idx) = false;
end