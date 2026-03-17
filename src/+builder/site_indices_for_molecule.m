function idx = site_indices_for_molecule(sys, uniqueMolID)
%SITE_INDICES_FOR_MOLECULE Return site indices belonging to one unique molecule.
%
% Input
%   sys          working system struct
%   uniqueMolID  scalar unique molecule ID
%
% Output
%   idx          column vector of site indices

if ~isfield(sys, 'site_mol_id') || isempty(sys.site_mol_id)
    error('sys.site_mol_id is missing or empty.');
end

idx = find(sys.site_mol_id == uniqueMolID);

end