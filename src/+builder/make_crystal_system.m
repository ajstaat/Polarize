function sys = make_crystal_system(crystal, model, opts)
%MAKE_CRYSTAL_SYSTEM Build the working system struct for one calculation.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end

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

    % Template/provenance metadata from the UNIT CELL
    sys.base_mol_id = [];
    if isfield(crystal, 'base_mol_id')
        sys.base_mol_id = crystal.base_mol_id;
    end

    sys.site_label = {};
    if isfield(crystal, 'site_label')
        sys.site_label = crystal.site_label;
    end

    sys.site_type = {};
    if isfield(crystal, 'site_type')
        sys.site_type = crystal.site_type;
    end

    sys.site_class = {};
    if isfield(crystal, 'site_class')
        sys.site_class = crystal.site_class;
    end

    if ~isempty(sys.base_mol_id), sys.base_mol_id = sys.base_mol_id(:); end
    if ~isempty(sys.site_label),  sys.site_label  = sys.site_label(:);  end
    if ~isempty(sys.site_type),   sys.site_type   = sys.site_type(:);   end
    if ~isempty(sys.site_class),  sys.site_class  = sys.site_class(:);  end

    % -----------------------------
    % Run-specific setup
    % -----------------------------
    sys.supercell_size = [1 1 1];
    if isfield(opts, 'supercell_size') && ~isempty(opts.supercell_size)
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

    bondScale = 1.20;
    if isfield(opts, 'bondScale') && ~isempty(opts.bondScale)
        bondScale = opts.bondScale;
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

    % New class-based typing
    sys.polarizable_classes = {};
    if isfield(model, 'polarizable_classes')
        sys.polarizable_classes = model.polarizable_classes;
    end

    sys.alpha_by_class = struct();
    if isfield(model, 'alpha_by_class')
        sys.alpha_by_class = model.alpha_by_class;
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
            error('Cannot build cart_coords because lattice is missing.');
        end
        sys.cart_coords = geom.frac_to_cart(sys.frac_coords, sys.lattice);
    end

    if isempty(sys.frac_coords) && ~isempty(sys.cart_coords)
        if isempty(sys.lattice)
            error('Cannot build frac_coords because lattice is missing.');
        end
        sys.frac_coords = geom.cart_to_frac(sys.cart_coords, sys.lattice);
    end

    validate_unit_cell_inputs(sys);

    nUnit = size(sys.cart_coords, 1);

    % -----------------------------
    % Unit-cell site properties
    % -----------------------------
    sys.site_pos = sys.cart_coords;

    sys.site_is_polarizable = false(nUnit, 1);
    sys.site_alpha = zeros(nUnit, 1);
    sys.site_charge = zeros(nUnit, 1);

    for i = 1:nUnit
        assigned = false;

        % Prefer class-based assignment if available
        if ~isempty(sys.site_class) && ~isempty(sys.polarizable_classes)
            thisClass = sys.site_class{i};

            if ismember(thisClass, sys.polarizable_classes)
                if ~isfield(sys.alpha_by_class, thisClass)
                    error('builder:make_crystal_system:MissingAlphaByClass', ...
                        'No alpha_by_class entry found for polarizable class "%s".', thisClass);
                end

                sys.site_is_polarizable(i) = true;
                sys.site_alpha(i) = sys.alpha_by_class.(thisClass);
                assigned = true;
            end
        end

        % Fallback to element-based assignment if class-based typing is absent
        if ~assigned
            thisType = sys.site_type{i};

            if ismember(thisType, sys.polarizable_types)
                if ~isfield(sys.alpha_by_type, thisType)
                    error('builder:make_crystal_system:MissingAlphaByType', ...
                        'No alpha_by_type entry found for polarizable site type "%s".', thisType);
                end

                sys.site_is_polarizable(i) = true;
                sys.site_alpha(i) = sys.alpha_by_type.(thisType);
            end
        end
    end

    sys.n_sites = nUnit;

    % -----------------------------
    % Build supercell blindly
    % -----------------------------
    unit = struct();
    unit.lattice = sys.lattice;
    unit.frac_coords = sys.frac_coords;
    unit.cart_coords = sys.cart_coords;

    unit.mol_id = sys.base_mol_id;   % provenance only
    unit.site_label = sys.site_label;
    unit.site_type = sys.site_type;
    unit.site_class = sys.site_class;

    unit.site_is_polarizable = sys.site_is_polarizable;
    unit.site_alpha = sys.site_alpha;
    unit.site_charge = sys.site_charge;

    super = geom.build_supercell(unit, sys.supercell_size);

    sys.unit_cell = unit;
    sys.supercell = super;

    sys.super_lattice = super.lattice;
    sys.site_pos = super.cart_coords;
    sys.site_frac = super.frac_coords;
    sys.site_label = super.site_label;
    sys.site_type = super.site_type;
    sys.site_class = super.site_class;
    sys.site_is_polarizable = super.site_is_polarizable;
    sys.site_alpha = super.site_alpha;
    sys.site_charge = super.site_charge;
    sys.cell_shift = super.cell_shift;
    sys.image_id = super.image_id;
    sys.unit_site_index = super.unit_site_index;
    sys.n_unit_sites = super.n_unit_sites;
    sys.n_cells = super.n_cells;
    sys.n_sites = super.n_sites;

    sys.base_mol_id = super.mol_id;

    % -----------------------------
    % Identify molecules on the built supercell itself
    % -----------------------------
    molinfo = builder.identify_supercell_molecules(super, bondScale, ...
        'Verbose', sys.verbose, ...
        'ProgressInterval', 250);

    sys.site_mol_id = molinfo.site_mol_id;
    sys.molecule_table = molinfo.component_table;
    sys.supercell_pbc_bond_graph = molinfo.A_pbc;

    % -----------------------------
    % Initialize active bookkeeping
    % -----------------------------
    sys.site_is_active = false(sys.n_sites, 1);
    sys.active_site_indices = {};

    % -----------------------------
    % Optional workflow operations
    % -----------------------------
    if ~isempty(sys.removed_molecules)
        sys = builder.remove_molecules(sys, sys.removed_molecules);
    end

    if ~isempty(sys.active_molecules)
        sys = builder.select_active_molecules(sys, sys.active_molecules);
    end

    sys.units = struct();
    sys.units.length = 'angstrom';
    sys.units.alpha  = 'angstrom^3';
    sys.units.charge = 'e';

end


function validate_unit_cell_inputs(sys)
    if isempty(sys.lattice) || ~isequal(size(sys.lattice), [3 3])
        error('builder:make_crystal_system:BadLattice', ...
            'A valid 3x3 lattice is required.');
    end

    if isempty(sys.frac_coords) || size(sys.frac_coords, 2) ~= 3
        error('builder:make_crystal_system:BadFracCoords', ...
            'frac_coords must be N x 3.');
    end

    if isempty(sys.cart_coords) || size(sys.cart_coords, 2) ~= 3
        error('builder:make_crystal_system:BadCartCoords', ...
            'cart_coords must be N x 3.');
    end

    if size(sys.frac_coords, 1) ~= size(sys.cart_coords, 1)
        error('builder:make_crystal_system:CoordCountMismatch', ...
            'frac_coords and cart_coords must contain the same number of sites.');
    end

    nSites = size(sys.cart_coords, 1);

    if isempty(sys.base_mol_id) || numel(sys.base_mol_id) ~= nSites
        error('builder:make_crystal_system:BadBaseMolId', ...
            'base_mol_id must have one entry per site.');
    end

    if isempty(sys.site_type) || numel(sys.site_type) ~= nSites
        error('builder:make_crystal_system:BadSiteType', ...
            'site_type must have one entry per site.');
    end

    if isempty(sys.site_label) || numel(sys.site_label) ~= nSites
        error('builder:make_crystal_system:BadSiteLabel', ...
            'site_label must have one entry per site.');
    end

    if isempty(sys.site_class) || numel(sys.site_class) ~= nSites
        error('builder:make_crystal_system:BadSiteClass', ...
            'site_class must have one entry per site.');
    end
end