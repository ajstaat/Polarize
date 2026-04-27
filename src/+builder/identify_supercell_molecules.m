function molinfo = identify_supercell_molecules(super, bondScale, varargin)
%IDENTIFY_SUPERCELL_MOLECULES Identify molecule components in a supercell.
%
% molinfo = builder.identify_supercell_molecules(super)
% molinfo = builder.identify_supercell_molecules(super, bondScale)
%
% This function intentionally uses a PBC bond graph on the replicated
% supercell. That reconstructs molecules that cross internal replicated-cell
% boundaries. It then checks each connected component in displayed Cartesian
% coordinates to decide whether that component is complete/selectable inside
% the displayed supercell box.
%
% Inputs
%   super      struct from geom.build_supercell
%   bondScale  covalent-radius bond scale, default 1.20
%
% Required super fields:
%   .frac_coords   N x 3 fractional coordinates in supercell basis
%   .cart_coords   N x 3 displayed Cartesian coordinates
%   .lattice       3 x 3 supercell lattice, rows are direct vectors
%   .site_type     N x 1 species labels
%   .units.length  'bohr' or 'angstrom' recommended
%
% Output
%   molinfo.site_mol_id      N x 1 molecule component ID
%   molinfo.component_table  struct table-like component metadata
%   molinfo.A_pbc            supercell PBC bond graph

if nargin < 2 || isempty(bondScale)
    bondScale = 1.20;
end

p = inputParser;
addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Method', 'auto', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});

verbose = p.Results.Verbose;
method = p.Results.Method;

validate_super(super);

nSites = size(super.cart_coords, 1);

if verbose
    fprintf('identify_supercell_molecules: building PBC graph on supercell...\n');
end

% io.build_pbc_bond_graph expects lattice distances compatible with
% io.bond_graph_tools covalent radii, which are in Angstrom.
toAng = local_length_to_angstrom_factor(super);

latticeAng = super.lattice * toAng;

[A_pbc, ~] = io.build_pbc_bond_graph( ...
    super.site_type, ...
    super.frac_coords, ...
    latticeAng, ...
    bondScale, ...
    'Verbose', verbose, ...
    'Method', method);

G = graph(A_pbc, 'upper');
siteMolID = conncomp(G).';

molIDs = unique(siteMolID, 'stable');
nMol = numel(molIDs);

component_table = struct();

component_table.molecule_id = molIDs(:);
component_table.unique_mol_id = molIDs(:);  % compatibility alias

component_table.site_indices = cell(nMol, 1);
component_table.n_sites = zeros(nMol, 1);
component_table.com = zeros(nMol, 3);

component_table.is_complete_in_display = false(nMol, 1);
component_table.n_display_fragments = zeros(nMol, 1);
component_table.largest_fragment_fraction = zeros(nMol, 1);

component_table.base_mol_id_mode = zeros(nMol, 1);
component_table.base_mol_id_unique_count = zeros(nMol, 1);

if isfield(super, 'cell_shift') && ~isempty(super.cell_shift)
    component_table.cell_shift_min = zeros(nMol, 3);
    component_table.cell_shift_max = zeros(nMol, 3);
end

sysLike = struct();
sysLike.site_pos = super.cart_coords;
sysLike.site_type = local_to_cell_column(super.site_type);
sysLike.site_mol_id = siteMolID;
sysLike.unique_mol_id = siteMolID;
sysLike.units = local_get_units(super);

for m = 1:nMol
    molID = molIDs(m);
    idx = find(siteMolID == molID);

    component_table.site_indices{m} = idx;
    component_table.n_sites(m) = numel(idx);
    component_table.com(m, :) = mean(super.cart_coords(idx, :), 1);

    displayInfo = builder.molecule_display_status(sysLike, molID, ...
        'BondScale', bondScale);

    component_table.is_complete_in_display(m) = displayInfo.is_complete_in_display;
    component_table.n_display_fragments(m) = displayInfo.n_fragments;
    component_table.largest_fragment_fraction(m) = displayInfo.largest_fragment_fraction;

    if isfield(super, 'mol_id') && ~isempty(super.mol_id)
        baseIDs = super.mol_id(idx);
        uBase = unique(baseIDs, 'stable');

        component_table.base_mol_id_unique_count(m) = numel(uBase);
        component_table.base_mol_id_mode(m) = mode(baseIDs);
    end

    if isfield(super, 'cell_shift') && ~isempty(super.cell_shift)
        component_table.cell_shift_min(m, :) = min(super.cell_shift(idx, :), [], 1);
        component_table.cell_shift_max(m, :) = max(super.cell_shift(idx, :), [], 1);
    end

    if verbose && (mod(m, 250) == 0 || m == nMol)
        fprintf('  checked display completeness for %d / %d molecules\n', m, nMol);
    end
end

molinfo = struct();
molinfo.site_mol_id = siteMolID;
molinfo.component_table = component_table;
molinfo.A_pbc = A_pbc;

if verbose
    fprintf('identify_supercell_molecules: n components = %d\n', nMol);
    fprintf('identify_supercell_molecules: complete in display = %d\n', ...
        nnz(component_table.is_complete_in_display));
end

end

% =========================================================================
% Local helpers
% =========================================================================

function validate_super(super)

required = {
    'lattice'
    'frac_coords'
    'cart_coords'
    'site_type'
};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(super, name) || isempty(super.(name))
        error('builder:identify_supercell_molecules:MissingField', ...
            'super.%s is required and missing/empty.', name);
    end
end

n = size(super.cart_coords, 1);

if size(super.cart_coords, 2) ~= 3 || ...
        size(super.frac_coords, 1) ~= n || ...
        size(super.frac_coords, 2) ~= 3
    error('builder:identify_supercell_molecules:BadCoordinates', ...
        'super.frac_coords and super.cart_coords must be N x 3.');
end

if ~isequal(size(super.lattice), [3 3])
    error('builder:identify_supercell_molecules:BadLattice', ...
        'super.lattice must be 3 x 3.');
end

if numel(super.site_type) ~= n
    error('builder:identify_supercell_molecules:BadSiteType', ...
        'super.site_type must have one entry per site.');
end

end

function units = local_get_units(super)

if isfield(super, 'units') && ~isempty(super.units)
    units = super.units;
else
    units = struct();
    units.length = 'bohr';
end

end

function factor = local_length_to_angstrom_factor(super)

BOHR2ANG = 1 / 1.8897259886;

units = local_get_units(super);

if ~isfield(units, 'length') || isempty(units.length)
    unit = 'bohr';
else
    unit = units.length;
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
        error('builder:identify_supercell_molecules:UnsupportedLengthUnit', ...
            'Unsupported super.units.length: %s', unit);
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
    error('builder:identify_supercell_molecules:BadTextMetadata', ...
        'Text metadata must be string/cellstr/cell.');
end

end