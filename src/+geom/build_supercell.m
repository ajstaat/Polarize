function super = build_supercell(unit, supercell_size)
%BUILD_SUPERCELL Replicate a unit cell into a supercell.
%
% Input
%   unit            struct with at least:
%                   .lattice      3x3
%                   .frac_coords  N x 3
%                   .cart_coords  N x 3
%
%                   optional:
%                   .mol_id
%                   .site_label
%                   .site_type
%                   .site_is_polarizable
%                   .site_alpha
%                   .site_charge
%
%   supercell_size  [nx ny nz]
%
% Output
%   super           struct with replicated data

if ~isfield(unit, 'lattice') || isempty(unit.lattice)
    error('unit.lattice is required.');
end

if ~isfield(unit, 'frac_coords') || isempty(unit.frac_coords)
    error('unit.frac_coords is required.');
end

if ~isfield(unit, 'cart_coords') || isempty(unit.cart_coords)
    error('unit.cart_coords is required.');
end

lattice = unit.lattice;
frac0 = unit.frac_coords;
cart0 = unit.cart_coords;

nUnit = size(frac0, 1);
if size(cart0, 1) ~= nUnit
    error('unit.frac_coords and unit.cart_coords must have the same number of sites.');
end

nx = supercell_size(1);
ny = supercell_size(2);
nz = supercell_size(3);

shifts = geom.translation_grid(supercell_size);
nCells = size(shifts, 1);
nTot = nUnit * nCells;

% Supercell lattice:
% rows are lattice vectors, so scale each row independently
superLattice = [
    nx * lattice(1, :);
    ny * lattice(2, :);
    nz * lattice(3, :)
];

% Preallocate
superFrac = zeros(nTot, 3);
superCart = zeros(nTot, 3);
cell_shift = zeros(nTot, 3);
image_id = zeros(nTot, 1);
unit_site_index = zeros(nTot, 1);

% Optional fields
hasMolID = isfield(unit, 'mol_id') && ~isempty(unit.mol_id);
hasLabel = isfield(unit, 'site_label') && ~isempty(unit.site_label);
hasType  = isfield(unit, 'site_type') && ~isempty(unit.site_type);
hasPol   = isfield(unit, 'site_is_polarizable') && ~isempty(unit.site_is_polarizable);
hasAlpha = isfield(unit, 'site_alpha') && ~isempty(unit.site_alpha);
hasQ     = isfield(unit, 'site_charge') && ~isempty(unit.site_charge);

if hasMolID
    superMolID = zeros(nTot, 1);
else
    superMolID = [];
end

if hasLabel
    superLabel = cell(nTot, 1);
else
    superLabel = {};
end

if hasType
    superType = cell(nTot, 1);
else
    superType = {};
end

if hasPol
    superPol = false(nTot, 1);
else
    superPol = [];
end

if hasAlpha
    superAlpha = zeros(nTot, size(unit.site_alpha, 2));
else
    superAlpha = [];
end

if hasQ
    superQ = zeros(nTot, 1);
else
    superQ = [];
end

% Replicate
row0 = 1;
for icell = 1:nCells
    shift = shifts(icell, :);

    rows = row0:(row0 + nUnit - 1);

    % Fractional coords in the supercell basis
    % If a unit cell point is at f and image shift is s, then in the larger
    % supercell fractional coordinate basis it is (f + s) ./ [nx ny nz].
    superFrac(rows, :) = (frac0 + shift) ./ supercell_size;

    % Cartesian coords in ordinary Cartesian space
    shiftCart = shift * lattice;
    superCart(rows, :) = cart0 + shiftCart;

    cell_shift(rows, :) = repmat(shift, nUnit, 1);
    image_id(rows) = icell;
    unit_site_index(rows) = (1:nUnit).';

    if hasMolID
        superMolID(rows) = unit.mol_id(:);
    end

    if hasLabel
        superLabel(rows) = unit.site_label(:);
    end

    if hasType
        superType(rows) = unit.site_type(:);
    end

    if hasPol
        superPol(rows) = unit.site_is_polarizable(:);
    end

    if hasAlpha
        superAlpha(rows, :) = unit.site_alpha;
    end

    if hasQ
        superQ(rows) = unit.site_charge(:);
    end

    row0 = row0 + nUnit;
end

super = struct();
super.supercell_size = supercell_size;

super.unit_lattice = lattice;
super.lattice = superLattice;

super.frac_coords = superFrac;
super.cart_coords = superCart;

super.cell_shift = cell_shift;           % integer image shift [ix iy iz]
super.image_id = image_id;               % image counter
super.unit_site_index = unit_site_index; % which unit-cell site this came from

if hasMolID
    super.mol_id = superMolID;
else
    super.mol_id = [];
end

if hasLabel
    super.site_label = superLabel;
else
    super.site_label = {};
end

if hasType
    super.site_type = superType;
else
    super.site_type = {};
end

if hasPol
    super.site_is_polarizable = superPol;
else
    super.site_is_polarizable = [];
end

if hasAlpha
    super.site_alpha = superAlpha;
else
    super.site_alpha = [];
end

if hasQ
    super.site_charge = superQ;
else
    super.site_charge = [];
end

super.n_unit_sites = nUnit;
super.n_cells = nCells;
super.n_sites = nTot;

end