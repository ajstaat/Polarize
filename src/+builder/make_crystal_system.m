function sys = make_crystal_system(crystal, model, opts)
%MAKE_CRYSTAL_SYSTEM Build a replicated, molecule-aware crystal system.
%
% sys = builder.make_crystal_system(crystal)
% sys = builder.make_crystal_system(crystal, model)
% sys = builder.make_crystal_system(crystal, model, opts)
% sys = builder.make_crystal_system(crystal, supercell_size)
%
% Inputs
%   crystal : unit-cell chemical template from io.import_contcar_as_crystal
%
%   model : optional struct with fields
%       .thole_a
%
%       .polarizable_classes
%       .alpha_by_class
%
%       .polarizable_types
%       .alpha_by_type
%
%       .alpha_units
%           'angstrom^3' or 'atomic_unit'
%           default = 'angstrom^3'
%
%   opts : optional struct with fields
%       .supercell_size   [nx ny nz], default [1 1 1]
%       .bondScale        scalar, default 1.20
%       .verbose          logical, default false
%
% Output
%   sys : replicated supercell system
%
% Internal unit convention:
%   sys.site_pos       bohr
%   sys.super_lattice  bohr
%   sys.site_alpha     atomic units
%   sys.site_charge    elementary charge
%
% Molecule identity:
%   Canonical molecule identity is assigned by
%   builder.identify_supercell_molecules(super, ...), which uses a
%   supercell PBC bond graph to reconstruct molecular connectivity, then
%   checks each component for displayed-coordinate completeness.

if nargin < 2 || isempty(model)
    model = struct();
end

if nargin < 3 || isempty(opts)
    opts = struct();
end

% Convenience form:
%   builder.make_crystal_system(crystal, [nx ny nz])
if isnumeric(model) && numel(model) == 3
    opts.supercell_size = reshape(model, 1, 3);
    model = struct();
end

supercellSize = local_get_opt(opts, 'supercell_size', [1 1 1]);
bondScale = local_get_opt(opts, 'bondScale', 1.20);
verbose = local_get_opt(opts, 'verbose', false);

supercellSize = reshape(supercellSize, 1, 3);

if any(supercellSize < 1) || any(mod(supercellSize, 1) ~= 0)
    error('builder:make_crystal_system:BadSupercellSize', ...
        'opts.supercell_size must be a 1x3 vector of positive integers.');
end

validate_crystal_template(crystal);

% -------------------------------------------------------------------------
% Unit-cell site metadata
% -------------------------------------------------------------------------

nUnit = crystal.nSites;

siteType = local_to_cell_column(crystal.site_type);
siteClass = local_to_cell_column(crystal.site_class);
siteLabel = local_to_cell_column(crystal.site_label);

unitIsPol = false(nUnit, 1);
unitAlpha = zeros(nUnit, 1);
unitCharge = zeros(nUnit, 1);

alphaUnits = local_get_field(model, 'alpha_units', 'angstrom^3');
alphaFactor = local_alpha_to_au_factor(alphaUnits);

polarizableClasses = local_to_cell_column_allow_empty( ...
    local_get_field(model, 'polarizable_classes', {}));

polarizableTypes = local_to_cell_column_allow_empty( ...
    local_get_field(model, 'polarizable_types', {}));

alphaByClass = local_get_field(model, 'alpha_by_class', struct());
alphaByType = local_get_field(model, 'alpha_by_type', struct());

for i = 1:nUnit
    assigned = false;

    if ~isempty(polarizableClasses)
        cls = siteClass{i};

        if ismember(cls, polarizableClasses)
            if ~isfield(alphaByClass, cls)
                error('builder:make_crystal_system:MissingAlphaByClass', ...
                    'model.alpha_by_class.%s is required.', cls);
            end

            unitIsPol(i) = true;
            unitAlpha(i) = alphaFactor * alphaByClass.(cls);
            assigned = true;
        end
    end

    if ~assigned && ~isempty(polarizableTypes)
        typ = siteType{i};

        if ismember(typ, polarizableTypes)
            if ~isfield(alphaByType, typ)
                error('builder:make_crystal_system:MissingAlphaByType', ...
                    'model.alpha_by_type.%s is required.', typ);
            end

            unitIsPol(i) = true;
            unitAlpha(i) = alphaFactor * alphaByType.(typ);
            assigned = true;
        end
    end
end

% -------------------------------------------------------------------------
% Build replicated supercell
% -------------------------------------------------------------------------

unit = struct();

unit.lattice = crystal.lattice;          % bohr
unit.frac_coords = crystal.frac_coords;  % dimensionless
unit.cart_coords = crystal.cart_coords;  % bohr

% This is provenance from the unit-cell crystal template. It is not the
% final molecule ID used for charge assignment / neighbor selection.
unit.mol_id = crystal.base_mol_id(:);

unit.site_type = siteType;
unit.site_class = siteClass;
unit.site_label = siteLabel;

unit.site_is_polarizable = unitIsPol;
unit.site_alpha = unitAlpha;             % atomic units
unit.site_charge = unitCharge;           % elementary charge

super = geom.build_supercell(unit, supercellSize);

% geom.build_supercell does not know physical units, so stamp them here for
% builder.identify_supercell_molecules.
super.units = struct();
super.units.length = 'bohr';
super.units.alpha = 'atomic_unit';
super.units.charge = 'elementary_charge';

% -------------------------------------------------------------------------
% Pack system fields before molecule identification
% -------------------------------------------------------------------------

sys = struct();

sys.comment = local_get_field(crystal, 'comment', '');
sys.source_file = local_get_field(crystal, 'source_file', '');

sys.units = struct();
sys.units.length = 'bohr';
sys.units.alpha = 'atomic_unit';
sys.units.charge = 'elementary_charge';

sys.thole_a = local_get_field(model, 'thole_a', []);

sys.cellpar = local_get_field(crystal, 'cellpar', []);
sys.unit_lattice = crystal.lattice;
sys.super_lattice = super.lattice;
sys.lattice = super.lattice;

sys.supercell_size = supercellSize;

sys.n_unit_sites = nUnit;
sys.n_base_mols = crystal.nBaseMols;
sys.n_cells = super.n_cells;
sys.n_sites = super.n_sites;

sys.site_pos = super.cart_coords;     % bohr
sys.site_frac = super.frac_coords;    % fractional in supercell basis

sys.cell_shift = super.cell_shift;
sys.image_id = super.image_id;
sys.unit_site_index = super.unit_site_index;

% Provenance fields from replicated unit-cell template.
sys.base_mol_id = super.mol_id(:);

sys.site_type = local_to_cell_column(super.site_type);
sys.site_class = local_to_cell_column(super.site_class);
sys.site_label = local_to_cell_column(super.site_label);

sys.site_is_polarizable = logical(super.site_is_polarizable(:));
sys.site_alpha = super.site_alpha(:);  % atomic units
sys.site_charge = super.site_charge(:);

% -------------------------------------------------------------------------
% Identify physical molecular components in the supercell
% -------------------------------------------------------------------------
% This is intentionally not [base_mol_id, cell_shift]. The supercell PBC
% graph reconstructs molecules that cross internal replicated-cell
% boundaries. The molecule table then records which components are complete
% in displayed Cartesian coordinates.

molinfo = builder.identify_supercell_molecules(super, bondScale, ...
    'Verbose', verbose);

sys.site_mol_id = molinfo.site_mol_id(:);
sys.unique_mol_id = sys.site_mol_id;  % compatibility alias

sys.molecule_table = molinfo.component_table;
sys.supercell_pbc_bond_graph = molinfo.A_pbc;

% Compatibility aliases for downstream code.
if isfield(sys.molecule_table, 'molecule_id') && ...
        ~isfield(sys.molecule_table, 'unique_mol_id')
    sys.molecule_table.unique_mol_id = sys.molecule_table.molecule_id;
elseif isfield(sys.molecule_table, 'unique_mol_id') && ...
        ~isfield(sys.molecule_table, 'molecule_id')
    sys.molecule_table.molecule_id = sys.molecule_table.unique_mol_id;
end

% Active/charged molecule bookkeeping filled in later by builder utilities.
sys.site_is_active = false(sys.n_sites, 1);
sys.active_site_indices = {};
sys.active_molecules = [];
sys.removed_molecules = [];

% Final consistency checks.
io.assert_atomic_units(sys);
validate_built_system(sys);

if verbose
    nMol = numel(sys.molecule_table.molecule_id);
    nComplete = nnz(sys.molecule_table.is_complete_in_display);

    fprintf('make_crystal_system: n unit sites        = %d\n', sys.n_unit_sites);
    fprintf('make_crystal_system: supercell size      = [%d %d %d]\n', supercellSize);
    fprintf('make_crystal_system: n supercell sites   = %d\n', sys.n_sites);
    fprintf('make_crystal_system: n molecules         = %d\n', nMol);
    fprintf('make_crystal_system: complete molecules  = %d\n', nComplete);
end

end

% =========================================================================
% Local validation
% =========================================================================

function validate_crystal_template(crystal)

required = {
    'lattice'
    'frac_coords'
    'cart_coords'
    'base_mol_id'
    'site_type'
    'site_class'
    'site_label'
    'units'
    'nSites'
    'nBaseMols'
};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(crystal, name) || isempty(crystal.(name))
        error('builder:make_crystal_system:MissingCrystalField', ...
            'crystal.%s is required and missing/empty.', name);
    end
end

if ~isequal(size(crystal.lattice), [3 3])
    error('builder:make_crystal_system:BadLattice', ...
        'crystal.lattice must be 3 x 3.');
end

if size(crystal.frac_coords, 2) ~= 3 || size(crystal.cart_coords, 2) ~= 3
    error('builder:make_crystal_system:BadCoordinates', ...
        'crystal.frac_coords and crystal.cart_coords must be N x 3.');
end

n = size(crystal.cart_coords, 1);

if size(crystal.frac_coords, 1) ~= n
    error('builder:make_crystal_system:CoordCountMismatch', ...
        'crystal.frac_coords and crystal.cart_coords must have same row count.');
end

if crystal.nSites ~= n
    error('builder:make_crystal_system:BadNSites', ...
        'crystal.nSites must match coordinate row count.');
end

if numel(crystal.base_mol_id) ~= n || ...
        numel(crystal.site_type) ~= n || ...
        numel(crystal.site_class) ~= n || ...
        numel(crystal.site_label) ~= n
    error('builder:make_crystal_system:BadMetadataLength', ...
        'Crystal site metadata must have one entry per site.');
end

if ~isfield(crystal.units, 'length') || ...
        ~strcmpi(crystal.units.length, 'bohr')
    error('builder:make_crystal_system:BadUnits', ...
        'crystal.units.length must be bohr.');
end

end

function validate_built_system(sys)

n = sys.n_sites;

required = {
    'site_pos'
    'site_frac'
    'site_type'
    'site_class'
    'site_label'
    'site_alpha'
    'site_charge'
    'site_is_polarizable'
    'site_mol_id'
    'molecule_table'
};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(sys, name) || isempty(sys.(name))
        error('builder:make_crystal_system:MissingBuiltField', ...
            'sys.%s is required and missing/empty.', name);
    end
end

if size(sys.site_pos, 1) ~= n || size(sys.site_pos, 2) ~= 3
    error('builder:make_crystal_system:BadSitePos', ...
        'sys.site_pos must be n_sites x 3.');
end

if size(sys.site_frac, 1) ~= n || size(sys.site_frac, 2) ~= 3
    error('builder:make_crystal_system:BadSiteFrac', ...
        'sys.site_frac must be n_sites x 3.');
end

if numel(sys.site_type) ~= n || ...
        numel(sys.site_class) ~= n || ...
        numel(sys.site_label) ~= n || ...
        numel(sys.site_alpha) ~= n || ...
        numel(sys.site_charge) ~= n || ...
        numel(sys.site_is_polarizable) ~= n || ...
        numel(sys.site_mol_id) ~= n
    error('builder:make_crystal_system:BadSiteFieldLength', ...
        'All site fields must have one entry per site.');
end

T = sys.molecule_table;

if ~isfield(T, 'molecule_id') || ...
        ~isfield(T, 'site_indices') || ...
        ~isfield(T, 'is_complete_in_display')
    error('builder:make_crystal_system:BadMoleculeTable', ...
        ['sys.molecule_table must contain molecule_id, site_indices, ' ...
         'and is_complete_in_display.']);
end

end

% =========================================================================
% Local helpers
% =========================================================================

function factor = local_alpha_to_au_factor(alphaUnits)

A3_PER_AU_ALPHA = 0.148184711;

if isstring(alphaUnits)
    alphaUnits = char(alphaUnits);
end

if ~ischar(alphaUnits)
    error('builder:make_crystal_system:BadAlphaUnits', ...
        'model.alpha_units must be a string scalar or character vector.');
end

alphaUnits = lower(strtrim(alphaUnits));

switch alphaUnits
    case {'angstrom^3', 'angstrom3', 'ang^3', 'ang3', 'a^3', 'aa^3'}
        factor = 1 / A3_PER_AU_ALPHA;

    case {'atomic_unit', 'atomic_units', 'au', 'a.u.', 'au_alpha', ...
            'bohr^3', 'bohr3'}
        factor = 1;

    otherwise
        error('builder:make_crystal_system:UnsupportedAlphaUnits', ...
            'Unsupported model.alpha_units: %s', alphaUnits);
end

end

function value = local_get_opt(s, name, defaultValue)

if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end

end

function value = local_get_field(s, name, defaultValue)

if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end

end

function c = local_to_cell_column(x)

if isempty(x)
    c = {};
elseif isstring(x)
    c = cellstr(x(:));
elseif iscellstr(x)
    c = x(:);
elseif iscell(x)
    c = x(:);
else
    error('builder:make_crystal_system:BadTextMetadata', ...
        'Text metadata must be string/cellstr/cell.');
end

end

function c = local_to_cell_column_allow_empty(x)

if isempty(x)
    c = {};
else
    c = local_to_cell_column(x);
end

end