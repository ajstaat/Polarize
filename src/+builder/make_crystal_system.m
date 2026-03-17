function sys = make_crystal_system(crystal, model, opts)
%MAKE_CRYSTAL_SYSTEM Build the working system struct for one calculation.
%
% This version:
%   - normalizes unit-cell geometry
%   - assigns site polarizabilities by type
%   - builds the replicated supercell working system
%   - assigns unique molecule IDs to each replicated molecule image
%   - optionally removes selected molecules
%   - optionally marks up to two active molecules

sys = struct();

% -----------------------------
% Copy basic crystal information
% -----------------------------
sys.cellpar = [];
if isfield(crystal, 'cellpar')
    sys.cellpar = crystal.cellpar;
end

sys.lattice = [];
if isfield(crystal, 'lattice')
    sys.lattice = crystal.lattice;
end

sys.frac_coords = [];
if isfield(crystal, 'frac_coords')
    sys.frac_coords = crystal.frac_coords;
end

sys.cart_coords = [];
if isfield(crystal, 'cart_coords')
    sys.cart_coords = crystal.cart_coords;
end

sys.mol_id = [];
if isfield(crystal, 'mol_id')
    sys.mol_id = crystal.mol_id;
end

sys.site_label = {};
if isfield(crystal, 'site_label')
    sys.site_label = crystal.site_label;
end

sys.site_type = {};
if isfield(crystal, 'site_type')
    sys.site_type = crystal.site_type;
end

% -----------------------------
% Run-specific setup
% -----------------------------
sys.supercell_size = [1 1 1];
if isfield(opts, 'supercell_size')
    sys.supercell_size = opts.supercell_size;
end

sys.active_molecules = [];
if isfield(opts, 'activeMolIDs')
    sys.active_molecules = opts.activeMolIDs;
end

sys.removed_molecules = [];
if isfield(opts, 'removeMolIDs')
    sys.removed_molecules = opts.removeMolIDs;
end

sys.verbose = false;
if isfield(opts, 'verbose')
    sys.verbose = opts.verbose;
end

% -----------------------------
% Model information
% -----------------------------
sys.thole_a = [];
if isfield(model, 'thole_a')
    sys.thole_a = model.thole_a;
end

sys.polarizable_types = {};
if isfield(model, 'polarizable_types')
    sys.polarizable_types = model.polarizable_types;
end

sys.alpha_by_type = struct();
if isfield(model, 'alpha_by_type')
    sys.alpha_by_type = model.alpha_by_type;
end

% -----------------------------
% Normalize unit-cell geometry
% -----------------------------
if isempty(sys.lattice)
    if ~isempty(sys.cellpar)
        sys.lattice = geom.lattice_vectors_from_cell(sys.cellpar);
    end
end

if isempty(sys.cart_coords) && ~isempty(sys.frac_coords)
    if isempty(sys.lattice)
        error(['Cannot build cart_coords from frac_coords because lattice is missing. ', ...
               'Provide crystal.lattice or crystal.cellpar.']);
    end
    sys.cart_coords = geom.frac_to_cart(sys.frac_coords, sys.lattice);
end

if isempty(sys.frac_coords) && ~isempty(sys.cart_coords)
    if isempty(sys.lattice)
        error(['Cannot build frac_coords from cart_coords because lattice is missing. ', ...
               'Provide crystal.lattice or crystal.cellpar.']);
    end
    sys.frac_coords = geom.cart_to_frac(sys.cart_coords, sys.lattice);
end

% -----------------------------
% Basic consistency checks
% -----------------------------
if ~isempty(sys.frac_coords) && size(sys.frac_coords, 2) ~= 3
    error('frac_coords must be N x 3.');
end

if ~isempty(sys.cart_coords) && size(sys.cart_coords, 2) ~= 3
    error('cart_coords must be N x 3.');
end

if ~isempty(sys.frac_coords) && ~isempty(sys.cart_coords)
    if size(sys.frac_coords, 1) ~= size(sys.cart_coords, 1)
        error('frac_coords and cart_coords must contain the same number of sites.');
    end
end

if ~isempty(sys.mol_id)
    nSites = size(sys.cart_coords, 1);
    if numel(sys.mol_id) ~= nSites
        error('mol_id must have one entry per site.');
    end
    sys.mol_id = sys.mol_id(:);
end

% -----------------------------
% Unit-cell site properties
% -----------------------------
nUnit = size(sys.cart_coords, 1);

sys.site_pos = sys.cart_coords;
sys.site_mol_id = sys.mol_id;

sys.site_is_polarizable = false(nUnit, 1);
sys.site_alpha = zeros(nUnit, 1);
sys.site_charge = zeros(nUnit, 1);

if ~isempty(sys.site_type)
    for i = 1:nUnit
        thisType = sys.site_type{i};

        if ismember(thisType, sys.polarizable_types)
            sys.site_is_polarizable(i) = true;

            if isfield(sys.alpha_by_type, thisType)
                sys.site_alpha(i) = sys.alpha_by_type.(thisType);
            else
                error('No alpha_by_type entry found for polarizable site type "%s".', thisType);
            end
        end
    end
end

sys.n_sites = nUnit;

% -----------------------------
% Build supercell working system
% -----------------------------
unit = struct();
unit.lattice = sys.lattice;
unit.frac_coords = sys.frac_coords;
unit.cart_coords = sys.cart_coords;
unit.mol_id = sys.mol_id;
unit.site_label = sys.site_label;
unit.site_type = sys.site_type;
unit.site_is_polarizable = sys.site_is_polarizable;
unit.site_alpha = sys.site_alpha;
unit.site_charge = sys.site_charge;

super = geom.build_supercell(unit, sys.supercell_size);

sys.unit_cell = unit;
sys.supercell = super;

sys.super_lattice = super.lattice;

sys.site_pos = super.cart_coords;
sys.site_frac = super.frac_coords;
sys.base_mol_id = super.mol_id;
sys.site_label = super.site_label;
sys.site_type = super.site_type;
sys.site_is_polarizable = super.site_is_polarizable;
sys.site_alpha = super.site_alpha;
sys.site_charge = super.site_charge;
sys.cell_shift = super.cell_shift;
sys.image_id = super.image_id;
sys.unit_site_index = super.unit_site_index;
sys.n_unit_sites = super.n_unit_sites;
sys.n_cells = super.n_cells;
sys.n_sites = super.n_sites;

[sys.unique_mol_id, sys.molecule_table] = ...
    builder.assign_unique_molecule_ids(sys.base_mol_id, sys.cell_shift);

sys.site_mol_id = sys.unique_mol_id;

% Initialize active bookkeeping
sys.site_is_active = false(sys.n_sites, 1);
sys.active_site_indices = {};

% -----------------------------
% Optional workflow operations
% -----------------------------
removeIDs = [];
if isfield(opts, 'removeMolIDs')
    removeIDs = opts.removeMolIDs;
end

activeIDs = [];
if isfield(opts, 'activeMolIDs')
    activeIDs = opts.activeMolIDs;
end

if ~isempty(removeIDs)
    sys = builder.remove_molecules(sys, removeIDs);
end

if ~isempty(activeIDs)
    sys = builder.select_active_molecules(sys, activeIDs);
end

end