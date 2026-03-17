function uniqueMolID = find_unique_molecule_id(sys, baseMolID, cellShift)
%FIND_UNIQUE_MOLECULE_ID Find unique molecule ID from baseMolID and cell shift.
%
% Input
%   sys        working system struct
%   baseMolID  scalar base molecule ID from the unit cell
%   cellShift  1x3 integer vector [ix iy iz]
%
% Output
%   uniqueMolID   scalar unique molecule ID
%
% Errors if no match or multiple matches are found.

if ~isfield(sys, 'molecule_table') || isempty(sys.molecule_table)
    error('sys.molecule_table is missing or empty.');
end

if numel(cellShift) ~= 3
    error('cellShift must be a 1x3 vector.');
end

cellShift = reshape(cellShift, 1, 3);

matchBase = (sys.molecule_table.base_mol_id == baseMolID);
matchShift = all(sys.molecule_table.cell_shift == cellShift, 2);

match = find(matchBase & matchShift);

if isempty(match)
    error('No molecule found for baseMolID = %d and cellShift = [%d %d %d].', ...
        baseMolID, cellShift(1), cellShift(2), cellShift(3));
end

if numel(match) > 1
    error('Multiple molecules found for baseMolID = %d and cellShift = [%d %d %d].', ...
        baseMolID, cellShift(1), cellShift(2), cellShift(3));
end

uniqueMolID = sys.molecule_table.unique_mol_id(match);

end