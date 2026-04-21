function molinfo = identify_supercell_molecules(super, bondScale, varargin)
%IDENTIFY_SUPERCELL_MOLECULES Identify molecules on a built supercell.
%
% molinfo = builder.identify_supercell_molecules(super)
% molinfo = builder.identify_supercell_molecules(super, bondScale)
% molinfo = builder.identify_supercell_molecules(super, bondScale, ...)
%
% Inputs
%   super      struct from geom.build_supercell, with at least:
%                .lattice
%                .frac_coords
%                .cart_coords
%                .site_type
%
%   bondScale  scalar, default 1.20
%
% Optional name-value inputs
%   'Verbose'           logical, default false
%   'ProgressInterval'  positive integer, default 2500
%
% Output
%   molinfo struct with fields:
%       .site_mol_id
%       .component_table
%       .A_pbc
%
% component_table fields
%   .molecule_id
%   .site_indices
%   .n_sites
%   .com
%   .is_complete_in_display
%   .n_display_fragments
%   .largest_fragment_fraction

    if nargin < 2 || isempty(bondScale)
        bondScale = 1.20;
    end

    p = inputParser;
    addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ProgressInterval', 2500, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    parse(p, varargin{:});
    opt = p.Results;

    validate_super(super);

    if opt.Verbose
        fprintf('identify_supercell_molecules: building PBC graph on supercell...\n');
    end

    [A_pbc, ~] = io.build_pbc_bond_graph( ...
        super.site_type, ...
        super.frac_coords, ...
        super.lattice, ...
        bondScale, ...
        'Verbose', opt.Verbose, ...
        'ProgressInterval', opt.ProgressInterval, ...
        'ReturnDcart', false);

    siteMolID = conncomp(graph(A_pbc)).';
    molIDs = unique(siteMolID, 'stable');
    nMol = numel(molIDs);

    if opt.Verbose
        fprintf('identify_supercell_molecules: found %d periodic molecules\n', nMol);
        fprintf('identify_supercell_molecules: checking display completeness...\n');
    end

    component_table = struct();
    component_table.molecule_id = molIDs(:);
    component_table.site_indices = cell(nMol, 1);
    component_table.n_sites = zeros(nMol, 1);
    component_table.com = zeros(nMol, 3);
    component_table.is_complete_in_display = false(nMol, 1);
    component_table.n_display_fragments = zeros(nMol, 1);
    component_table.largest_fragment_fraction = zeros(nMol, 1);

    sysLike = struct();
    sysLike.site_pos = super.cart_coords;
    sysLike.site_type = super.site_type;
    sysLike.site_mol_id = siteMolID;

    for m = 1:nMol
        thisMolID = molIDs(m);
        idx = find(siteMolID == thisMolID);

        displayInfo = builder.molecule_display_status(sysLike, thisMolID, ...
            'BondScale', bondScale);

        component_table.site_indices{m} = idx;
        component_table.n_sites(m) = numel(idx);
        component_table.com(m, :) = mean(super.cart_coords(idx, :), 1);
        component_table.is_complete_in_display(m) = displayInfo.is_complete_in_display;
        component_table.n_display_fragments(m) = displayInfo.n_fragments;
        component_table.largest_fragment_fraction(m) = displayInfo.largest_fragment_fraction;

        if opt.Verbose
            if mod(m, opt.ProgressInterval) == 0 || m == nMol
                fprintf('  checked %d / %d molecules\n', m, nMol);
            end
        end
    end

    if opt.Verbose
        nComplete = nnz(component_table.is_complete_in_display);
        fragMin = min(component_table.n_display_fragments);
        fragMax = max(component_table.n_display_fragments);
        fprintf('identify_supercell_molecules: complete molecules = %d / %d\n', ...
            nComplete, nMol);
        fprintf('identify_supercell_molecules: display fragment count range = [%d, %d]\n', ...
            fragMin, fragMax);
    end

    molinfo = struct();
    molinfo.site_mol_id = siteMolID;
    molinfo.component_table = component_table;
    molinfo.A_pbc = A_pbc;
end


function validate_super(super)
    required = {'lattice', 'frac_coords', 'cart_coords', 'site_type'};
    for k = 1:numel(required)
        name = required{k};
        if ~isfield(super, name) || isempty(super.(name))
            error('builder:identify_supercell_molecules:MissingField', ...
                'super.%s is required and missing/empty.', name);
        end
    end

    if ~isequal(size(super.lattice), [3 3])
        error('builder:identify_supercell_molecules:BadLattice', ...
            'super.lattice must be 3 x 3.');
    end

    if size(super.frac_coords, 2) ~= 3 || size(super.cart_coords, 2) ~= 3
        error('builder:identify_supercell_molecules:BadCoords', ...
            'super.frac_coords and super.cart_coords must be N x 3.');
    end

    if size(super.frac_coords, 1) ~= size(super.cart_coords, 1)
        error('builder:identify_supercell_molecules:CoordCountMismatch', ...
            'super.frac_coords and super.cart_coords must have the same number of rows.');
    end

    if numel(super.site_type) ~= size(super.cart_coords, 1)
        error('builder:identify_supercell_molecules:BadSiteType', ...
            'super.site_type must have one entry per site.');
    end
end