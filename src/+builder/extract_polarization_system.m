function polsys = extract_polarization_system(sys, params)
%EXTRACT_POLARIZATION_SYSTEM Build canonical solver-facing polarization system.
%
% This function converts a rich crystal/system struct into a canonical
% polarization-system struct ("polsys") suitable for both nonperiodic and
% periodic operator assembly.
%
% Responsibilities:
%   - keep site-level polarization data
%   - preserve lattice/supercell information needed for periodic assembly
%   - preserve useful site metadata
%   - normalize units to internal atomic units
%   - attach provenance/meta information
%
% Required fields on input sys:
%   .n_sites
%   .site_pos
%   .site_alpha
%   .site_charge
%   .site_is_polarizable
%   .thole_a
%
% Optional but preserved if present:
%   .site_mol_id
%   .site_label
%   .site_type
%   .lattice
%   .super_lattice
%   .cellpar
%   .supercell_size
%   .active_molecules
%   .units
%
% Inputs
%   sys     rich system struct
%   params  calculation params (used mainly for periodic mode metadata)
%
% Output
%   polsys  canonical polarization system in atomic units

    if nargin < 2 || isempty(params)
        params = struct();
    end

    % ---------------------------------------------------------------------
    % Validate required fields
    % ---------------------------------------------------------------------
    required = { ...
        'n_sites', ...
        'site_pos', ...
        'site_alpha', ...
        'site_charge', ...
        'site_is_polarizable', ...
        'thole_a'};

    for k = 1:numel(required)
        if ~isfield(sys, required{k}) || isempty(sys.(required{k}))
            error('sys.%s is required and missing/empty.', required{k});
        end
    end

    nSites = sys.n_sites;

    if size(sys.site_pos, 1) ~= nSites || size(sys.site_pos, 2) ~= 3
        error('sys.site_pos must be n_sites x 3.');
    end
    if numel(sys.site_alpha) ~= nSites
        error('sys.site_alpha must have length n_sites.');
    end
    if numel(sys.site_charge) ~= nSites
        error('sys.site_charge must have length n_sites.');
    end
    if numel(sys.site_is_polarizable) ~= nSites
        error('sys.site_is_polarizable must have length n_sites.');
    end

    % ---------------------------------------------------------------------
    % Start canonical struct by copying only relevant fields
    % ---------------------------------------------------------------------
    polsys = struct();

    % Core site-level polarization data
    polsys.n_sites = nSites;
    polsys.site_pos = sys.site_pos;
    polsys.site_alpha = sys.site_alpha(:);
    polsys.site_charge = sys.site_charge(:);
    polsys.site_is_polarizable = logical(sys.site_is_polarizable(:));
    polsys.thole_a = sys.thole_a;

    % Optional site metadata
    if isfield(sys, 'site_mol_id') && ~isempty(sys.site_mol_id)
        polsys.site_mol_id = sys.site_mol_id(:);
    else
        polsys.site_mol_id = (1:nSites).';
    end

    if isfield(sys, 'site_label') && ~isempty(sys.site_label)
        polsys.site_label = sys.site_label;
    else
        polsys.site_label = arrayfun(@(k) sprintf('site_%d', k), 1:nSites, 'UniformOutput', false).';
    end

    if isfield(sys, 'site_type') && ~isempty(sys.site_type)
        polsys.site_type = sys.site_type;
    else
        polsys.site_type = repmat({'X'}, nSites, 1);
    end

    % ---------------------------------------------------------------------
    % Preserve lattice / supercell geometry for periodic calculations
    % ---------------------------------------------------------------------
    polsys.lattice = [];
    if isfield(sys, 'lattice') && ~isempty(sys.lattice)
        polsys.lattice = sys.lattice;
    end

    polsys.super_lattice = [];
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        polsys.super_lattice = sys.super_lattice;
    end

    polsys.cellpar = [];
    if isfield(sys, 'cellpar') && ~isempty(sys.cellpar)
        polsys.cellpar = sys.cellpar;
    end

    polsys.supercell_size = [];
    if isfield(sys, 'supercell_size') && ~isempty(sys.supercell_size)
        polsys.supercell_size = sys.supercell_size;
    end

    % Optional workflow metadata
    if isfield(sys, 'active_molecules') && ~isempty(sys.active_molecules)
        polsys.active_molecules = sys.active_molecules(:);
    else
        polsys.active_molecules = [];
    end

    % ---------------------------------------------------------------------
    % Periodic metadata
    % ---------------------------------------------------------------------
    polsys.periodic = struct();
    polsys.periodic.mode = '';
    if isfield(params, 'ewald') && isfield(params.ewald, 'mode') && ~isempty(params.ewald.mode)
        polsys.periodic.mode = char(string(params.ewald.mode));
    end

    polsys.periodic.is_periodic = strcmpi(polsys.periodic.mode, 'periodic_triclinic');

    H = [];
    if ~isempty(polsys.super_lattice)
        H = polsys.super_lattice;
    elseif ~isempty(polsys.lattice)
        H = polsys.lattice;
    end

    polsys.periodic.real_space_lattice = H;
    polsys.periodic.recip_lattice = [];
    polsys.periodic.volume = [];

    if ~isempty(H)
        if ~isequal(size(H), [3, 3])
            error('lattice/super_lattice must be 3x3 when present.');
        end

        % Reciprocal lattice convention:
        %   B = 2*pi * inv(H)'
        polsys.periodic.volume = abs(det(H));
        polsys.periodic.recip_lattice = 2*pi * inv(H).';
    end

    % ---------------------------------------------------------------------
    % Units: default legacy assumption is Angstrom / Angstrom^3 / e
    % ---------------------------------------------------------------------
    if ~isfield(polsys, 'units') || isempty(polsys)
        % no-op, just to keep structure clear
    end

    if isfield(sys, 'units') && ~isempty(sys.units)
        polsys.units = sys.units;
    else
        polsys.units = struct();
        polsys.units.length = 'angstrom';
        polsys.units.alpha = 'angstrom^3';
        polsys.units.charge = 'e';
    end

    % ---------------------------------------------------------------------
    % Provenance / metadata
    % ---------------------------------------------------------------------
    polsys.meta = struct();
    polsys.meta.source = 'builder.extract_polarization_system';
    polsys.meta.original_n_sites = nSites;
    polsys.meta.n_polarizable_sites = nnz(polsys.site_is_polarizable);

    if isfield(sys, 'meta') && ~isempty(sys.meta)
        polsys.meta.upstream_meta = sys.meta;
    end

    % ---------------------------------------------------------------------
    % Normalize to atomic units once, here
    % ---------------------------------------------------------------------
    polsys = io.convert_to_atomic_units(polsys);

    % ---------------------------------------------------------------------
    % Add a simple source mask for nonzero charges
    % ---------------------------------------------------------------------
    polsys.source_mask = abs(polsys.site_charge) > 0;
    polsys.target_mask = true(nSites, 1);
end