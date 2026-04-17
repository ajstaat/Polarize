function crystal = import_contcar_as_crystal(filename, varargin)
%IMPORT_CONTCAR_AS_CRYSTAL Read a VASP POSCAR/CONTCAR-style file and
% return a crystal template struct for Polarize.
%
% This function builds the UNIT-CELL CHEMICAL TEMPLATE:
%   - wrapped unit-cell coordinates
%   - base molecule IDs found in the unit cell
%   - local site labels within each base molecule
%   - site types
%   - simple graph-based site classes for polarizability typing
%
% It does NOT define final molecule identities for the working supercell.
% Those are identified later on the built supercell itself.

    p = inputParser;
    addRequired(p, 'filename', @(x) ischar(x) || isstring(x));
    addParameter(p, 'BondScale', 1.20, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'WrapFractional', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'SortMolecules', false, @(x) islogical(x) && isscalar(x));
    parse(p, filename, varargin{:});
    opts = p.Results;

    S = io.read_vasp_structure(filename);

    if opts.WrapFractional
        S.frac = S.frac - floor(S.frac);
        S.cart = S.frac * S.lattice;
    end

    % Unit-cell PBC bond graph
    [A_unit, ~] = io.build_pbc_bond_graph(S.species, S.frac, S.lattice, opts.BondScale);

    % Identify molecules in the UNIT CELL for template/provenance purposes
    molecules = io.unwrap_all_contcar_molecules(S, ...
        'BondScale', opts.BondScale, ...
        'SortMolecules', opts.SortMolecules);

    nSites = size(S.frac, 1);
    base_mol_id = zeros(nSites, 1);

    for k = 1:numel(molecules)
        idx = molecules{k}.indices(:);
        base_mol_id(idx) = k;
    end

    if any(base_mol_id == 0)
        error('io:import_contcar_as_crystal:UnassignedSites', ...
            'Failed to assign base molecule IDs to all sites.');
    end

    % New graph-based classes
    site_class = io.infer_simple_site_classes(S.species, A_unit);

    crystal = struct();
    crystal.cellpar     = [];
    crystal.lattice     = S.lattice;
    crystal.frac_coords = S.frac;
    crystal.cart_coords = S.cart;

    crystal.base_mol_id = base_mol_id;
    crystal.site_type   = S.species(:);
    crystal.site_class  = site_class(:);
    crystal.site_label  = io.make_molecule_local_site_labels(S.species, base_mol_id);

    validate_crystal_import(crystal);
end

function validate_crystal_import(crystal)
    n = size(crystal.cart_coords, 1);

    if ~isequal(size(crystal.lattice), [3 3])
        error('io:import_contcar_as_crystal:BadLattice', ...
            'crystal.lattice must be 3x3.');
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
end