function info = molecule_display_status(sys, uniqueMolID, varargin)
%MOLECULE_DISPLAY_STATUS Check whether one molecule component is displayed complete.
%
% info = builder.molecule_display_status(sys, uniqueMolID)
% info = builder.molecule_display_status(sys, uniqueMolID, 'BondScale', 1.20)
%
% This checks connectivity using displayed Cartesian coordinates only.
% If a molecule component is split across the displayed supercell boundary,
% it will have multiple displayed fragments.
%
% Distance-unit convention:
%   sys.site_pos may be in bohr or Angstrom. io.bond_graph_tools uses
%   covalent radii in Angstrom, so distances are converted to Angstrom
%   before applying the bond test.

p = inputParser;
addRequired(p, 'sys', @isstruct);
addRequired(p, 'uniqueMolID', @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'BondScale', 1.20, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
parse(p, sys, uniqueMolID, varargin{:});

bondScale = p.Results.BondScale;

validate_sys(sys);

idx = builder.site_indices_for_molecule(sys, uniqueMolID);

if isempty(idx)
    error('builder:molecule_display_status:UnknownMolID', ...
        'uniqueMolID %g was not found.', uniqueMolID);
end

X = sys.site_pos(idx, :);
species = local_to_cell_column(sys.site_type(idx));

toAng = local_length_to_angstrom_factor(sys);

n = size(X, 1);

if n == 1
    fragID = 1;
else
    A = false(n, n);

    for i = 1:(n - 1)
        si = species{i};

        for j = (i + 1):n
            sj = species{j};

            d = norm(X(j, :) - X(i, :)) * toAng;

            if io.bond_graph_tools.is_bonded(si, sj, d, bondScale)
                A(i, j) = true;
                A(j, i) = true;
            end
        end
    end

    fragID = conncomp(graph(A)).';
end

fragIDs = unique(fragID, 'stable');
nFrag = numel(fragIDs);

fragSizes = zeros(nFrag, 1);

for k = 1:nFrag
    fragSizes(k) = nnz(fragID == fragIDs(k));
end

info = struct();

info.unique_mol_id = uniqueMolID;
info.molecule_id = uniqueMolID;

info.site_indices = idx(:);
info.n_sites = n;

info.fragment_ids = fragID(:);
info.fragment_sizes = fragSizes(:);
info.n_fragments = nFrag;

info.is_complete_in_display = (nFrag == 1);
info.largest_fragment_fraction = max(fragSizes) / n;
info.com = mean(X, 1);

end

% =========================================================================
% Local helpers
% =========================================================================

function validate_sys(sys)

required = {'site_pos', 'site_type'};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(sys, name) || isempty(sys.(name))
        error('builder:molecule_display_status:MissingField', ...
            'sys.%s is required and missing/empty.', name);
    end
end

if size(sys.site_pos, 2) ~= 3
    error('builder:molecule_display_status:BadSitePos', ...
        'sys.site_pos must be N x 3.');
end

if numel(sys.site_type) ~= size(sys.site_pos, 1)
    error('builder:molecule_display_status:BadSiteType', ...
        'sys.site_type must have one entry per site.');
end

if ~(isfield(sys, 'site_mol_id') && ~isempty(sys.site_mol_id)) && ...
        ~(isfield(sys, 'unique_mol_id') && ~isempty(sys.unique_mol_id))
    error('builder:molecule_display_status:MissingMoleculeIDs', ...
        'sys.site_mol_id or sys.unique_mol_id is required.');
end

end

function factor = local_length_to_angstrom_factor(sys)

BOHR2ANG = 1 / 1.8897259886;

if ~isfield(sys, 'units') || isempty(sys.units) || ...
        ~isfield(sys.units, 'length') || isempty(sys.units.length)
    % Rebuilt systems should have units.length = 'bohr'. This default is
    % chosen for safety in the current refactor branch.
    unit = 'bohr';
else
    unit = sys.units.length;
end

if isstring(unit)
    unit = char(unit);
end

unit = lower(strtrim(unit));

switch unit
    case {'angstrom', 'angstroms', 'ang', 'a'}
        factor = 1.0;

    case {'bohr', 'bohrs', 'a0', 'au_length'}
        factor = BOHR2ANG;

    otherwise
        error('builder:molecule_display_status:UnsupportedLengthUnit', ...
            'Unsupported sys.units.length: %s', unit);
end

end

function c = local_to_cell_column(x)

if isstring(x)
    c = cellstr(x(:));
elseif iscellstr(x)
    c = x(:);
elseif iscell(x)
    c = x(:);
else
    error('builder:molecule_display_status:BadTextMetadata', ...
        'Text metadata must be string/cellstr/cell.');
end

end