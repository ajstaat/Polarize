function super = build_supercell(unit, supercell_size)
%BUILD_SUPERCELL Replicate a unit cell into a supercell.
%
% Project lattice convention:
%   lattice rows are direct lattice vectors
%   cart = frac * lattice
%
% Required input fields:
%   unit.lattice       3 x 3 direct lattice, rows are lattice vectors
%   unit.frac_coords   N x 3 fractional coordinates in the unit cell
%   unit.cart_coords   N x 3 Cartesian coordinates in the unit cell
%
% Optional site fields, replicated if present:
%   unit.mol_id
%   unit.site_label
%   unit.site_type
%   unit.site_class
%   unit.site_is_polarizable
%   unit.site_alpha
%   unit.site_charge
%
% Input:
%   supercell_size     [nx ny nz], positive integer replication counts
%
% Output:
%   super              struct with replicated coordinates, lattice, image
%                      metadata, and replicated optional site fields.

% -------------------------------------------------------------------------
% Validate required fields
% -------------------------------------------------------------------------

if ~isfield(unit, 'lattice') || isempty(unit.lattice)
    error('geom:build_supercell:MissingLattice', ...
        'unit.lattice is required.');
end

if ~isfield(unit, 'frac_coords') || isempty(unit.frac_coords)
    error('geom:build_supercell:MissingFracCoords', ...
        'unit.frac_coords is required.');
end

if ~isfield(unit, 'cart_coords') || isempty(unit.cart_coords)
    error('geom:build_supercell:MissingCartCoords', ...
        'unit.cart_coords is required.');
end

lattice = unit.lattice;
frac0 = unit.frac_coords;
cart0 = unit.cart_coords;

if ~isnumeric(lattice) || ~isequal(size(lattice), [3 3])
    error('geom:build_supercell:BadLattice', ...
        'unit.lattice must be a numeric 3x3 matrix.');
end

if ~isnumeric(frac0) || size(frac0, 2) ~= 3
    error('geom:build_supercell:BadFracCoords', ...
        'unit.frac_coords must be a numeric N x 3 array.');
end

if ~isnumeric(cart0) || size(cart0, 2) ~= 3
    error('geom:build_supercell:BadCartCoords', ...
        'unit.cart_coords must be a numeric N x 3 array.');
end

nUnit = size(frac0, 1);

if size(cart0, 1) ~= nUnit
    error('geom:build_supercell:CoordCountMismatch', ...
        'unit.frac_coords and unit.cart_coords must have the same number of sites.');
end

if ~isnumeric(supercell_size) || numel(supercell_size) ~= 3 || ...
        any(supercell_size < 1) || any(mod(supercell_size, 1) ~= 0)
    error('geom:build_supercell:BadSupercellSize', ...
        'supercell_size must be a 1x3 vector of positive integers.');
end

supercell_size = reshape(supercell_size, 1, 3);

% -------------------------------------------------------------------------
% Validate optional fields
% -------------------------------------------------------------------------

hasMolID = local_has_field(unit, 'mol_id');
hasLabel = local_has_field(unit, 'site_label');
hasType  = local_has_field(unit, 'site_type');
hasClass = local_has_field(unit, 'site_class');
hasPol   = local_has_field(unit, 'site_is_polarizable');
hasAlpha = local_has_field(unit, 'site_alpha');
hasQ     = local_has_field(unit, 'site_charge');

if hasMolID
    local_assert_vector_length(unit.mol_id, nUnit, 'unit.mol_id');
end

if hasLabel
    local_assert_text_vector_length(unit.site_label, nUnit, 'unit.site_label');
end

if hasType
    local_assert_text_vector_length(unit.site_type, nUnit, 'unit.site_type');
end

if hasClass
    local_assert_text_vector_length(unit.site_class, nUnit, 'unit.site_class');
end

if hasPol
    local_assert_vector_length(unit.site_is_polarizable, nUnit, ...
        'unit.site_is_polarizable');
end

if hasAlpha
    if ~isnumeric(unit.site_alpha) || size(unit.site_alpha, 1) ~= nUnit
        error('geom:build_supercell:BadSiteAlpha', ...
            'unit.site_alpha must be a numeric array with one row per unit-cell site.');
    end
end

if hasQ
    local_assert_vector_length(unit.site_charge, nUnit, 'unit.site_charge');
end

% -------------------------------------------------------------------------
% Construct supercell lattice and preallocate required fields
% -------------------------------------------------------------------------

nx = supercell_size(1);
ny = supercell_size(2);
nz = supercell_size(3);

shifts = geom.translation_grid(supercell_size);
nCells = size(shifts, 1);
nTot = nUnit * nCells;

% Rows are direct lattice vectors.
superLattice = [
    nx * lattice(1, :)
    ny * lattice(2, :)
    nz * lattice(3, :)
];

superFrac = zeros(nTot, 3);
superCart = zeros(nTot, 3);

cell_shift = zeros(nTot, 3);
image_id = zeros(nTot, 1);
unit_site_index = zeros(nTot, 1);

% -------------------------------------------------------------------------
% Preallocate optional fields, preserving text container type where sensible
% -------------------------------------------------------------------------

if hasMolID
    superMolID = zeros(nTot, 1);
else
    superMolID = [];
end

if hasLabel
    superLabel = local_prealloc_text_like(unit.site_label, nTot);
else
    superLabel = [];
end

if hasType
    superType = local_prealloc_text_like(unit.site_type, nTot);
else
    superType = [];
end

if hasClass
    superClass = local_prealloc_text_like(unit.site_class, nTot);
else
    superClass = [];
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

% -------------------------------------------------------------------------
% Replicate coordinates and metadata
% -------------------------------------------------------------------------

row0 = 1;

for icell = 1:nCells
    shift = shifts(icell, :);
    rows = row0:(row0 + nUnit - 1);

    % Fractional coordinates in the supercell basis.
    superFrac(rows, :) = (frac0 + shift) ./ supercell_size;

    % Cartesian coordinates in ordinary Cartesian space.
    shiftCart = shift * lattice;
    superCart(rows, :) = cart0 + shiftCart;

    cell_shift(rows, :) = repmat(shift, nUnit, 1);
    image_id(rows) = icell;
    unit_site_index(rows) = (1:nUnit).';

    if hasMolID
        superMolID(rows) = unit.mol_id(:);
    end

    if hasLabel
        superLabel(rows) = local_cast_text_like(unit.site_label(:), superLabel);
    end

    if hasType
        superType(rows) = local_cast_text_like(unit.site_type(:), superType);
    end

    if hasClass
        superClass(rows) = local_cast_text_like(unit.site_class(:), superClass);
    end

    if hasPol
        superPol(rows) = logical(unit.site_is_polarizable(:));
    end

    if hasAlpha
        superAlpha(rows, :) = unit.site_alpha;
    end

    if hasQ
        superQ(rows) = unit.site_charge(:);
    end

    row0 = row0 + nUnit;
end

% -------------------------------------------------------------------------
% Pack output
% -------------------------------------------------------------------------

super = struct();

super.supercell_size = supercell_size;

super.unit_lattice = lattice;
super.lattice = superLattice;

super.frac_coords = superFrac;
super.cart_coords = superCart;

super.cell_shift = cell_shift;           % integer image shift [ix iy iz]
super.image_id = image_id;               % image counter
super.unit_site_index = unit_site_index; % source unit-cell site index

if hasMolID
    super.mol_id = superMolID;
else
    super.mol_id = [];
end

if hasLabel
    super.site_label = superLabel;
else
    super.site_label = [];
end

if hasType
    super.site_type = superType;
else
    super.site_type = [];
end

if hasClass
    super.site_class = superClass;
else
    super.site_class = [];
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

% =========================================================================
% Local helpers
% =========================================================================

function tf = local_has_field(s, name)
tf = isfield(s, name) && ~isempty(s.(name));
end

function local_assert_vector_length(x, n, fieldName)
if numel(x) ~= n
    error('geom:build_supercell:BadOptionalFieldLength', ...
        '%s must have one entry per unit-cell site.', fieldName);
end
end

function local_assert_text_vector_length(x, n, fieldName)
if ~(isstring(x) || iscellstr(x) || ischar(x))
    error('geom:build_supercell:BadTextField', ...
        '%s must be a string array, cellstr, or char array.', fieldName);
end

if ischar(x)
    if size(x, 1) ~= n
        error('geom:build_supercell:BadTextFieldLength', ...
            '%s char array must have one row per unit-cell site.', fieldName);
    end
else
    if numel(x) ~= n
        error('geom:build_supercell:BadTextFieldLength', ...
            '%s must have one entry per unit-cell site.', fieldName);
    end
end
end

function out = local_prealloc_text_like(x, n)
if isstring(x)
    out = strings(n, 1);
elseif iscellstr(x)
    out = cell(n, 1);
elseif ischar(x)
    % Preserve legacy char-array style as cellstr in the supercell. This is
    % easier to assign safely when labels have variable widths.
    out = cell(n, 1);
else
    error('geom:build_supercell:BadTextField', ...
        'Text metadata must be a string array, cellstr, or char array.');
end
end

function y = local_cast_text_like(x, target)
if isstring(target)
    y = string(x);
elseif iscell(target)
    if isstring(x)
        y = cellstr(x);
    elseif ischar(x)
        y = cellstr(x);
    elseif iscellstr(x)
        y = x;
    else
        error('geom:build_supercell:BadTextCast', ...
            'Cannot cast text metadata to cellstr.');
    end
else
    error('geom:build_supercell:BadTextTarget', ...
        'Unsupported text metadata target type.');
end

y = y(:);
end