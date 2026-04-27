function molIDs = complete_molecule_ids(sys, varargin)
%COMPLETE_MOLECULE_IDS Return molecule IDs complete in the displayed box.
%
% molIDs = builder.complete_molecule_ids(sys)
% molIDs = builder.complete_molecule_ids(sys, 'Verbose', true)
%
% This is the main selection gate for downstream reference/neighbor picking.
% A molecule ID is returned only if sys.molecule_table marks it complete in
% displayed Cartesian coordinates.

p = inputParser;
addRequired(p, 'sys', @isstruct);
addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
parse(p, sys, varargin{:});

verbose = p.Results.Verbose;

validate_sys(sys);

T = sys.molecule_table;

if isfield(T, 'molecule_id') && ~isempty(T.molecule_id)
    ids = T.molecule_id(:);
elseif isfield(T, 'unique_mol_id') && ~isempty(T.unique_mol_id)
    ids = T.unique_mol_id(:);
else
    error('builder:complete_molecule_ids:MissingMoleculeID', ...
        'sys.molecule_table must contain molecule_id or unique_mol_id.');
end

mask = logical(T.is_complete_in_display(:));

if numel(mask) ~= numel(ids)
    error('builder:complete_molecule_ids:SizeMismatch', ...
        'molecule_id and is_complete_in_display must have the same length.');
end

molIDs = ids(mask);
molIDs = molIDs(:);

if verbose
    fprintf('Complete molecule summary:\n');
    fprintf('  total molecules found        = %d\n', numel(ids));
    fprintf('  complete in displayed box    = %d\n', numel(molIDs));

    if isempty(molIDs)
        fprintf('  complete molecule IDs        = []\n');
    else
        fprintf('  complete molecule IDs        = %s\n', mat2str(molIDs(:).'));
    end
end

end

function validate_sys(sys)

if ~isfield(sys, 'molecule_table') || isempty(sys.molecule_table)
    error('builder:complete_molecule_ids:MissingMoleculeTable', ...
        'sys.molecule_table is required and missing/empty.');
end

T = sys.molecule_table;

hasID = (isfield(T, 'molecule_id') && ~isempty(T.molecule_id)) || ...
        (isfield(T, 'unique_mol_id') && ~isempty(T.unique_mol_id));

if ~hasID
    error('builder:complete_molecule_ids:BadMoleculeTable', ...
        'sys.molecule_table must contain molecule_id or unique_mol_id.');
end

if ~isfield(T, 'is_complete_in_display') || isempty(T.is_complete_in_display)
    error('builder:complete_molecule_ids:BadMoleculeTable', ...
        'sys.molecule_table.is_complete_in_display is required.');
end

end