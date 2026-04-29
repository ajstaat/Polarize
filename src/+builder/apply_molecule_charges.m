function sys = apply_molecule_charges(sys, molIDs, varargin)
%APPLY_MOLECULE_CHARGES Assign charges to selected molecules.
%
% sys = builder.apply_molecule_charges(sys, molIDs, Name, Value)
%
% Common uniform-charge usage:
%
%   sys = builder.apply_molecule_charges(sys, [refID nbrID], ...
%       'Mode', 'uniform', ...
%       'TotalCharges', [+1 -1], ...
%       'SetActive', true, ...
%       'DisablePolarizabilityOnCharged', true, ...
%       'ZeroExistingCharges', true);
%
% Inputs
%   sys     system struct from builder.make_crystal_system
%   molIDs  molecule IDs to charge
%
% Name-value options
%   'Mode'
%       'uniform' or 'file'
%       default: 'uniform'
%
%   'TotalCharges'
%       vector of total molecular charges, one per molID.
%       For uniform mode, each molecule's total charge is distributed
%       equally over that molecule's sites.
%       default: ones(size(molIDs))
%
%   'SetActive'
%       if true, call builder.select_active_molecules(sys, molIDs)
%       default: true
%
%   'DisablePolarizabilityOnCharged'
%       if true, charged molecule sites are removed from the polarizable set
%       default: true
%
%   'ZeroExistingCharges'
%       if true, zero all sys.site_charge before assigning new charges
%       default: false
%
%   'RequireComplete'
%       if true, all molIDs must be complete in displayed coordinates
%       default: true
%
%   'ChargeFile'
%       file used in Mode='file'. Expected columns are:
%           label charge
%       or
%           x y z charge
%       depending on later workflow needs.
%       This mode is retained for compatibility but is intentionally kept
%       strict in this refactor branch.
%
%   'Template'
%       optional template struct for Mode='file' with fields:
%           .site_pos    N x 3
%           .site_type   N x 1
%           .site_charge N x 1
%
%   'DistanceTol'
%       matching tolerance for Mode='file', in same length units as
%       template.site_pos and sys.site_pos
%       default: 1e-3
%
% Units
%   sys.site_charge is in elementary charge.
%   sys.site_alpha is in atomic units and is not rescaled here.

p = inputParser;
addRequired(p, 'sys', @isstruct);
addRequired(p, 'molIDs', @(x) isnumeric(x) && isvector(x) && ~isempty(x));
addParameter(p, 'Mode', 'uniform', @(x) ischar(x) || isstring(x));
addParameter(p, 'TotalCharges', [], @(x) isnumeric(x) && isvector(x));
addParameter(p, 'SetActive', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'DisablePolarizabilityOnCharged', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ZeroExistingCharges', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'RequireComplete', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ChargeFile', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Template', struct(), @isstruct);
addParameter(p, 'DistanceTol', 0.2, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'ReferenceAxis', [0 0 1], @(x) isnumeric(x) && numel(x) == 3 && norm(x) > 0);
addParameter(p, 'PrimaryAxis', [1 0 0], @(x) isnumeric(x) && numel(x) == 3 && norm(x) > 0);
addParameter(p, 'AmbiguityTol', 1e-10, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
parse(p, sys, molIDs, varargin{:});

opt = p.Results;

molIDs = molIDs(:).';
mode = lower(char(string(opt.Mode)));

validate_sys(sys);

if isempty(opt.TotalCharges)
    totalCharges = ones(size(molIDs));
else
    totalCharges = reshape(opt.TotalCharges, 1, []);
end

if numel(totalCharges) ~= numel(molIDs)
    error('builder:apply_molecule_charges:ChargeCountMismatch', ...
        'TotalCharges must contain one value per molecule ID.');
end

if opt.RequireComplete
    assert_molecules_complete(sys, molIDs);
end

if opt.ZeroExistingCharges
    sys.site_charge(:) = 0;
end

switch mode
    case 'uniform'
        sys = assign_uniform_charges(sys, molIDs, totalCharges, opt.Verbose);

    case 'file'
        sys = assign_file_charges(sys, molIDs, totalCharges, opt);

    otherwise
        error('builder:apply_molecule_charges:BadMode', ...
            'Unsupported charge assignment mode: %s', mode);
end

if opt.DisablePolarizabilityOnCharged
    for k = 1:numel(molIDs)
        idx = builder.site_indices_for_molecule(sys, molIDs(k));

        % Remove charged molecule sites from the induced-dipole active
        % space, but preserve their physical site_alpha values.
        %
        % prepare_scf_problem uses site_is_polarizable to decide which
        % sites participate as induced dipoles, so zeroing site_alpha is
        % not needed for SCF exclusion.
        %
        % Keeping site_alpha is important for Thole-damped external fields:
        % charge -> polarizable-site damping uses both source and target
        % alpha values. If charged source-site alpha is zeroed here, then
        % fieldParams.use_thole_damping=true cannot actually apply the
        % short-range Thole correction for charged source sites.
        sys.site_is_polarizable(idx) = false;
    end
end

if opt.SetActive
    sys = builder.select_active_molecules(sys, molIDs);
end

sys.charged_molecules = molIDs(:);
sys.charged_molecule_total_charges = totalCharges(:);

if isfield(sys, 'units')
    sys.units.charge = 'elementary_charge';
end

end

% =========================================================================
% Assignment modes
% =========================================================================

function sys = assign_uniform_charges(sys, molIDs, totalCharges, verbose)

for k = 1:numel(molIDs)
    molID = molIDs(k);
    idx = builder.site_indices_for_molecule(sys, molID);

    if isempty(idx)
        error('builder:apply_molecule_charges:UnknownMolecule', ...
            'Molecule ID %d was not found.', molID);
    end

    qTotal = totalCharges(k);
    qSite = qTotal / numel(idx);

    sys.site_charge(idx) = qSite;

    if verbose
        fprintf('Assigned uniform charge to molecule %d:\n', molID);
        fprintf('  n sites      = %d\n', numel(idx));
        fprintf('  total charge = %+ .6f e\n', qTotal);
        fprintf('  site charge  = %+ .6f e\n', qSite);
    end
end

end

function sys = assign_file_charges(sys, molIDs, totalCharges, opt)
%ASSIGN_FILE_CHARGES Compatibility path for file/template charge mapping.
%
% This path is intentionally strict. The first supported refactor workflow is
% uniform molecular charging. File-based mapping is kept here so the API has
% an obvious landing spot, but robust tests for it belong in Builder-3b.

template = opt.Template;

if isempty(fieldnames(template))
    if strlength(string(opt.ChargeFile)) == 0
        error('builder:apply_molecule_charges:MissingTemplate', ...
            ['Mode=''file'' requires either Template or ChargeFile. ' ...
             'File-mode tests will be added in Builder-3b.']);
    end

    template = local_read_charge_template(opt.ChargeFile);
end

required = {'site_pos', 'site_type', 'site_charge'};
for r = 1:numel(required)
    name = required{r};

    if ~isfield(template, name) || isempty(template.(name))
        error('builder:apply_molecule_charges:BadTemplate', ...
            'Template must contain nonempty field "%s".', name);
    end
end

template.site_type = local_to_cell_column(template.site_type);
template.site_charge = template.site_charge(:);

if size(template.site_pos, 1) ~= numel(template.site_type) || ...
        numel(template.site_charge) ~= numel(template.site_type)
    error('builder:apply_molecule_charges:BadTemplateSize', ...
        'Template site_pos, site_type, and site_charge lengths must match.');
end

for k = 1:numel(molIDs)
    molID = molIDs(k);
    idx = builder.site_indices_for_molecule(sys, molID);

    target = struct();
    target.site_pos = sys.site_pos(idx, :);
    target.site_type = local_to_cell_column(sys.site_type(idx));

    map = builder.match_molecule_atoms_by_frame(template, target, ...
        'DistanceTol', opt.DistanceTol, ...
        'ReferenceAxis', opt.ReferenceAxis, ...
        'PrimaryAxis', opt.PrimaryAxis, ...
        'AmbiguityTol', opt.AmbiguityTol);
    
    % q must be in target-site order. target_to_template(j) gives the template
    % atom corresponding to target atom j.
    q = template.site_charge(map.target_to_template);

    % If TotalCharges was supplied, rescale template charges to requested
    % total molecular charge.
    qSum = sum(q);
    qTarget = totalCharges(k);

    if abs(qSum) > eps
        q = q * (qTarget / qSum);
    elseif abs(qTarget) > eps
        error('builder:apply_molecule_charges:CannotRescaleZeroTemplateCharge', ...
            'Template charges sum to zero but requested total charge is nonzero.');
    end

    sys.site_charge(idx) = q(:);
end

end

% =========================================================================
% Validation / helpers
% =========================================================================

function validate_sys(sys)

required = {
    'site_pos'
    'site_mol_id'
    'site_type'
    'site_charge'
    'site_is_polarizable'
    'site_alpha'
    'molecule_table'
};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(sys, name) || isempty(sys.(name))
        error('builder:apply_molecule_charges:MissingField', ...
            'sys.%s is required and missing/empty.', name);
    end
end

n = size(sys.site_pos, 1);

if numel(sys.site_mol_id) ~= n || ...
        numel(sys.site_type) ~= n || ...
        numel(sys.site_charge) ~= n || ...
        numel(sys.site_is_polarizable) ~= n || ...
        numel(sys.site_alpha) ~= n
    error('builder:apply_molecule_charges:BadSiteFieldLength', ...
        'Site fields must have one entry per site.');
end

T = sys.molecule_table;

if ~isfield(T, 'molecule_id') || isempty(T.molecule_id)
    error('builder:apply_molecule_charges:BadMoleculeTable', ...
        'sys.molecule_table.molecule_id is required.');
end

end

function assert_molecules_complete(sys, molIDs)

T = sys.molecule_table;

if ~isfield(T, 'is_complete_in_display') || isempty(T.is_complete_in_display)
    error('builder:apply_molecule_charges:MissingCompletenessFlag', ...
        'sys.molecule_table.is_complete_in_display is required.');
end

for k = 1:numel(molIDs)
    row = find(T.molecule_id == molIDs(k), 1, 'first');

    if isempty(row)
        error('builder:apply_molecule_charges:UnknownMolecule', ...
            'Molecule ID %d was not found.', molIDs(k));
    end

    if ~T.is_complete_in_display(row)
        error('builder:apply_molecule_charges:IncompleteMolecule', ...
            'Molecule ID %d is not complete in the displayed supercell.', molIDs(k));
    end
end

end

function template = local_read_charge_template(filename)
%LOCAL_READ_CHARGE_TEMPLATE Small compatibility reader.
%
% Expected numeric format:
%   x y z charge
%
% Element labels are not inferable from a pure numeric file, so this reader
% is deliberately conservative. Rich file-mode support belongs in Builder-3b.

filename = char(string(filename));

if ~isfile(filename)
    error('builder:apply_molecule_charges:MissingChargeFile', ...
        'Charge file not found: %s', filename);
end

M = readmatrix(filename);

if size(M, 2) < 4
    error('builder:apply_molecule_charges:BadChargeFile', ...
        'ChargeFile numeric format must contain at least four columns: x y z charge.');
end

template = struct();
template.site_pos = M(:, 1:3);
template.site_charge = M(:, 4);

error('builder:apply_molecule_charges:ChargeFileNeedsTypes', ...
    ['ChargeFile numeric reader found coordinates/charges, but site_type labels ' ...
     'are required for frame-based atom matching. Use Template for now.']);

end

function c = local_to_cell_column(x)

if isstring(x)
    c = cellstr(x(:));
elseif iscellstr(x)
    c = x(:);
elseif iscell(x)
    c = x(:);
else
    error('builder:apply_molecule_charges:BadTextMetadata', ...
        'Text metadata must be string/cellstr/cell.');
end

end