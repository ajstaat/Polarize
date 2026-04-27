function polsys = extract_polarization_system(sys, params)
%EXTRACT_POLARIZATION_SYSTEM Extract solver-facing polarization system.
%
% polsys = builder.extract_polarization_system(sys)
% polsys = builder.extract_polarization_system(sys, params)
%
% The builder system is already expected to be in atomic units:
%   site_pos      bohr
%   site_alpha    atomic units
%   site_charge   elementary charge
%
% This function keeps the active/charged molecules in the site list. The
% polarizable mask determines which sites participate as induced dipoles.
%
% params optional fields:
%   .mode
%       'periodic' or 'nonperiodic'
%
%   .ewald.mode
%       compatibility alias for params.mode
%
%   .include_inactive
%       default true
%
%   .active_only
%       if true, keep only sys.site_is_active sites
%       default false
%
% Output polsys fields:
%   site_pos
%   site_alpha
%   site_charge
%   site_is_polarizable
%   site_type
%   site_class
%   site_label
%   site_mol_id
%   site_is_active
%   units
%   thole_a
%   is_periodic
%   lattice / super_lattice

if nargin < 2 || isempty(params)
    params = struct();
end

validate_sys(sys);

mode = resolve_mode(params);
activeOnly = local_get_field(params, 'active_only', false);

if activeOnly
    if ~isfield(sys, 'site_is_active') || isempty(sys.site_is_active)
        error('builder:extract_polarization_system:MissingActiveMask', ...
            'active_only=true requires sys.site_is_active.');
    end

    keepMask = logical(sys.site_is_active(:));
else
    keepMask = true(size(sys.site_pos, 1), 1);
end

polsys = struct();

siteFields = {
    'site_pos'
    'site_alpha'
    'site_charge'
    'site_is_polarizable'
    'site_type'
    'site_class'
    'site_label'
    'site_mol_id'
    'unique_mol_id'
    'base_mol_id'
    'cell_shift'
    'image_id'
    'unit_site_index'
    'site_is_active'
};

for k = 1:numel(siteFields)
    name = siteFields{k};

    if isfield(sys, name) && ~isempty(sys.(name))
        polsys.(name) = local_slice_site_field(sys.(name), keepMask);
    end
end

polsys.n_sites = size(polsys.site_pos, 1);

polsys.units = sys.units;
polsys.units.length = 'bohr';
polsys.units.alpha = 'atomic_unit';
polsys.units.charge = 'elementary_charge';

polsys.thole_a = local_get_field(sys, 'thole_a', []);

if isempty(polsys.thole_a)
    warning('builder:extract_polarization_system:MissingTholeA', ...
        'sys.thole_a is missing/empty.');
end

if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
    polsys.super_lattice = sys.super_lattice;
    polsys.lattice = sys.super_lattice;
elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
    polsys.lattice = sys.lattice;
    polsys.super_lattice = sys.lattice;
end

if isfield(sys, 'supercell_size')
    polsys.supercell_size = sys.supercell_size;
end

switch mode
    case {'periodic', 'periodic_triclinic'}
        polsys.is_periodic = true;
        polsys.periodic_mode = 'periodic';

    case {'nonperiodic', 'finite', 'cluster'}
        polsys.is_periodic = false;
        polsys.periodic_mode = 'nonperiodic';

    otherwise
        error('builder:extract_polarization_system:BadMode', ...
            'Unsupported extraction mode: %s', mode);
end

% Keep molecule table if it exists, but do not try to reindex it for
% active_only extraction. The solver should use site masks, not molecule
% table rows.
if isfield(sys, 'molecule_table')
    polsys.molecule_table = sys.molecule_table;
end

io.assert_atomic_units(polsys);
validate_polsys(polsys);

end

% =========================================================================
% Validation / helpers
% =========================================================================

function validate_sys(sys)

required = {
    'site_pos'
    'site_alpha'
    'site_charge'
    'site_is_polarizable'
    'units'
};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(sys, name) || isempty(sys.(name))
        error('builder:extract_polarization_system:MissingField', ...
            'sys.%s is required and missing/empty.', name);
    end
end

n = size(sys.site_pos, 1);

if size(sys.site_pos, 2) ~= 3
    error('builder:extract_polarization_system:BadSitePos', ...
        'sys.site_pos must be N x 3.');
end

if numel(sys.site_alpha) ~= n || ...
        numel(sys.site_charge) ~= n || ...
        numel(sys.site_is_polarizable) ~= n
    error('builder:extract_polarization_system:BadSiteFieldLength', ...
        'site_alpha, site_charge, and site_is_polarizable must match site_pos.');
end

io.assert_atomic_units(sys);

end

function validate_polsys(polsys)

n = size(polsys.site_pos, 1);

if numel(polsys.site_alpha) ~= n || ...
        numel(polsys.site_charge) ~= n || ...
        numel(polsys.site_is_polarizable) ~= n
    error('builder:extract_polarization_system:BadPolsysLength', ...
        'Extracted polsys site fields have inconsistent lengths.');
end

if ~islogical(polsys.site_is_polarizable)
    polsys.site_is_polarizable = logical(polsys.site_is_polarizable);
end

end

function mode = resolve_mode(params)

mode = '';

if isfield(params, 'mode') && ~isempty(params.mode)
    mode = char(string(params.mode));
elseif isfield(params, 'ewald') && isstruct(params.ewald) && ...
        isfield(params.ewald, 'mode') && ~isempty(params.ewald.mode)
    mode = char(string(params.ewald.mode));
end

if isempty(mode)
    mode = 'periodic';
end

mode = lower(strtrim(mode));

end

function value = local_get_field(s, name, defaultValue)

if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end

end

function y = local_slice_site_field(x, keepMask)

if isvector(x) && numel(x) == numel(keepMask)
    y = x(keepMask);
elseif size(x, 1) == numel(keepMask)
    y = x(keepMask, :);
else
    y = x;
end

end