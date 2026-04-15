function crystal = import_contcar_as_crystal(filename, varargin)
%IMPORT_CONTCAR_AS_CRYSTAL Read a VASP POSCAR/CONTCAR-style file and
% return a Polarize-compatible crystal struct.
%
% Required output fields:
%   crystal.lattice
%   crystal.frac_coords
%   crystal.cart_coords
%   crystal.mol_id
%   crystal.site_label
%   crystal.site_type
%
% Optional name-value inputs:
%   'BondScale'      : scalar, default 1.20
%   'WrapFractional' : logical, default true
%   'SortMolecules'  : logical, default false

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

    [A, ~] = io.build_vasp_bond_graph(S, opts.BondScale);
    comp = conncomp(graph(A)).';

    if opts.SortMolecules
        comp = relabel_molecules_by_com(S, comp);
    end

    crystal = struct();
    crystal.cellpar     = [];
    crystal.lattice     = S.lattice;
    crystal.frac_coords = S.frac;
    crystal.cart_coords = S.cart;
    crystal.mol_id      = comp(:);
    crystal.site_type   = S.species(:);
    crystal.site_label  = io.make_molecule_local_site_labels(S.species, comp);

    validate_crystal_import(crystal);
end

function mol_id = relabel_molecules_by_com(S, mol_id_in)
    mols = unique(mol_id_in(:));
    coms = zeros(numel(mols), 3);

    for k = 1:numel(mols)
        idx = find(mol_id_in == mols(k));
        Mol = io.unwrap_vasp_molecule(S, idx);
        coms(k,:) = mean(Mol.cart, 1);
    end

    [~, order] = sortrows(coms, [3 2 1]);

    mol_id = zeros(size(mol_id_in));
    for newID = 1:numel(order)
        oldID = mols(order(newID));
        mol_id(mol_id_in == oldID) = newID;
    end
end

function validate_crystal_import(crystal)
    n = size(crystal.cart_coords, 1);

    if ~isequal(size(crystal.lattice), [3 3])
        error('crystal.lattice must be 3x3.');
    end
    if size(crystal.frac_coords, 2) ~= 3 || size(crystal.cart_coords, 2) ~= 3
        error('crystal.frac_coords and crystal.cart_coords must be N x 3.');
    end
    if size(crystal.frac_coords, 1) ~= n
        error('frac/cart coordinate count mismatch.');
    end
    if numel(crystal.mol_id) ~= n
        error('crystal.mol_id must have one entry per site.');
    end
    if numel(crystal.site_type) ~= n || numel(crystal.site_label) ~= n
        error('site_type and site_label must have one entry per site.');
    end
end