function crystal = import_contcar_as_crystal(filename, varargin)
%IMPORT_CONTCAR_AS_CRYSTAL Import a VASP POSCAR/CONTCAR as a crystal template.
%
% crystal = io.import_contcar_as_crystal(filename)
% crystal = io.import_contcar_as_crystal(filename, Name, Value)
%
% This function builds the UNIT-CELL CHEMICAL TEMPLATE:
%   - unit-cell lattice
%   - wrapped unit-cell coordinates
%   - base molecule IDs found in the unit cell
%   - local site labels within each base molecule
%   - simple graph-based site classes for polarizability typing
%
% It does NOT define final molecule identities for the working supercell.
% Those are identified later on the built supercell itself.
%
% Output unit convention:
%   crystal.lattice      bohr
%   crystal.cart_coords  bohr
%   crystal.frac_coords  dimensionless
%   crystal.units.length = 'bohr'

p = inputParser;
addRequired(p, 'filename', @(x) ischar(x) || isstring(x));
addParameter(p, 'BondScale', 1.20, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'WrapFractional', true, ...
    @(x) islogical(x) && isscalar(x));
addParameter(p, 'SortMolecules', false, ...
    @(x) islogical(x) && isscalar(x));
parse(p, filename, varargin{:});

opts = p.Results;

% -------------------------------------------------------------------------
% Read raw VASP structure
% -------------------------------------------------------------------------

S = io.read_vasp_structure(filename);

if opts.WrapFractional
    S.frac = S.frac - floor(S.frac);
    S.cart = S.frac * S.lattice;
end

nSites = size(S.frac, 1);

% -------------------------------------------------------------------------
% Unit-cell PBC bond graph
% -------------------------------------------------------------------------

[A_unit, ~] = io.build_pbc_bond_graph( ...
    S.species, S.frac, S.lattice, opts.BondScale);

% -------------------------------------------------------------------------
% Identify unit-cell connected components / base molecules
% -------------------------------------------------------------------------

molecules = io.unwrap_all_contcar_molecules(S, ...
    'BondScale', opts.BondScale, ...
    'BondGraph', A_unit, ...
    'SortMolecules', opts.SortMolecules);

base_mol_id = zeros(nSites, 1);

for k = 1:numel(molecules)
    idx = molecules{k}.indices(:);
    base_mol_id(idx) = k;
end

if any(base_mol_id == 0)
    error('io:import_contcar_as_crystal:UnassignedSites', ...
        'Failed to assign base molecule IDs to all sites.');
end

% -------------------------------------------------------------------------
% Metadata
% -------------------------------------------------------------------------

site_class = io.infer_simple_site_classes(S.species, A_unit);
site_label = io.make_molecule_local_site_labels(S.species, base_mol_id);

% -------------------------------------------------------------------------
% Convert lattice/cartesian coordinates to internal atomic units
% -------------------------------------------------------------------------

unitSys = struct();
unitSys.site_pos = S.cart;
unitSys.units.length = 'angstrom';
unitSys.units.alpha = 'atomic_unit';
unitSys.units.charge = 'elementary_charge';

unitSys = io.convert_to_atomic_units(unitSys);

latticeSys = struct();
latticeSys.site_pos = S.lattice;
latticeSys.units.length = 'angstrom';
latticeSys.units.alpha = 'atomic_unit';
latticeSys.units.charge = 'elementary_charge';

latticeSys = io.convert_to_atomic_units(latticeSys);

% -------------------------------------------------------------------------
% Pack crystal template
% -------------------------------------------------------------------------

crystal = struct();

crystal.comment = S.comment;
crystal.source_file = char(string(filename));

crystal.cellpar = [];
crystal.lattice = latticeSys.site_pos;
crystal.frac_coords = S.frac;
crystal.cart_coords = unitSys.site_pos;

crystal.base_mol_id = base_mol_id;

crystal.site_type = S.species(:);
crystal.site_class = site_class(:);
crystal.site_label = site_label(:);

crystal.units = struct();
crystal.units.length = 'bohr';
crystal.units.alpha = 'atomic_unit';
crystal.units.charge = 'elementary_charge';

crystal.nSites = nSites;
crystal.nBaseMols = numel(molecules);

validate_crystal_import(crystal);

end

function validate_crystal_import(crystal)

n = size(crystal.cart_coords, 1);

if ~isequal(size(crystal.lattice), [3 3])
    error('io:import_contcar_as_crystal:BadLattice', ...
        'crystal.lattice must be 3 x 3.');
end

if size(crystal.frac_coords, 2) ~= 3 || size(crystal.cart_coords, 2) ~= 3
    error('io:import_contcar_as_crystal:BadCoordinateShape', ...
        'crystal.frac_coords and crystal.cart_coords must be N x 3.');
end

if size(crystal.frac_coords, 1) ~= n
    error('io:import_contcar_as_crystal:CoordCountMismatch', ...
        'frac/cart coordinate count mismatch.');
end

if numel(crystal.base_mol_id) ~= n
    error('io:import_contcar_as_crystal:BadBaseMolIdSize', ...
        'crystal.base_mol_id must have one entry per site.');
end

if numel(crystal.site_type) ~= n || ...
        numel(crystal.site_class) ~= n || ...
        numel(crystal.site_label) ~= n
    error('io:import_contcar_as_crystal:BadSiteMetadataSize', ...
        'site_type, site_class, and site_label must have one entry per site.');
end

if ~isfield(crystal, 'units') || ...
        ~isfield(crystal.units, 'length') || ...
        ~strcmp(crystal.units.length, 'bohr')
    error('io:import_contcar_as_crystal:BadUnits', ...
        'crystal.units.length must be bohr.');
end

if ~isfield(crystal, 'nSites') || crystal.nSites ~= n
    error('io:import_contcar_as_crystal:BadNSites', ...
        'crystal.nSites must match the number of coordinate rows.');
end

if ~isfield(crystal, 'nBaseMols') || crystal.nBaseMols < 1
    error('io:import_contcar_as_crystal:BadNBaseMols', ...
        'crystal.nBaseMols must be positive.');
end

end