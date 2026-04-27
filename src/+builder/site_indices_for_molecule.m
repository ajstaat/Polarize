function idx = site_indices_for_molecule(sys, moleculeID)
%SITE_INDICES_FOR_MOLECULE Return site indices for one molecule component.
%
% idx = builder.site_indices_for_molecule(sys, moleculeID)
%
% Uses sys.site_mol_id as the canonical site-to-molecule assignment. This
% field is produced by builder.identify_supercell_molecules.

if ~isstruct(sys)
    error('builder:site_indices_for_molecule:BadSys', ...
        'sys must be a struct.');
end

if ~isscalar(moleculeID) || ~isnumeric(moleculeID) || ~isfinite(moleculeID)
    error('builder:site_indices_for_molecule:BadMoleculeID', ...
        'moleculeID must be a finite numeric scalar.');
end

if isfield(sys, 'site_mol_id') && ~isempty(sys.site_mol_id)
    molID = sys.site_mol_id(:);
elseif isfield(sys, 'unique_mol_id') && ~isempty(sys.unique_mol_id)
    molID = sys.unique_mol_id(:);
else
    error('builder:site_indices_for_molecule:MissingMoleculeIDs', ...
        'sys.site_mol_id or sys.unique_mol_id is required.');
end

idx = find(molID == moleculeID);
idx = idx(:);

end