function sys = convert_to_atomic_units(sys)
%CONVERT_TO_ATOMIC_UNITS Convert supported system fields to atomic units.
%
% Internal Polarize unit convention:
%   positions        : bohr
%   polarizability   : atomic units
%   charges          : elementary charge
%
% Expected fields, when present:
%   sys.site_pos       N x 3 positions
%   sys.site_alpha     N x 1, N x 3, or N x 3 x 3 polarizabilities
%   sys.site_charge    N x 1 charges
%
% Optional input metadata:
%   sys.units.length   'angstrom' or 'bohr'
%   sys.units.alpha    'angstrom^3', 'angstrom3', 'A^3', or 'atomic_unit'
%   sys.units.charge   'elementary_charge'
%
% If sys.units is missing, legacy defaults are assumed:
%   length = 'angstrom'
%   alpha  = 'angstrom^3'
%   charge = 'elementary_charge'
%
% Output:
%   sys with converted fields and stamped unit metadata.

% -------------------------------------------------------------------------
% Conversion constants
% -------------------------------------------------------------------------

ANG2BOHR = 1.8897259886;

% 1 atomic unit of polarizability = 0.148184711 Angstrom^3.
A3_PER_AU_ALPHA = 0.148184711;

% -------------------------------------------------------------------------
% Read or assign input unit metadata
% -------------------------------------------------------------------------

if ~isfield(sys, 'units') || isempty(sys.units)
    unitsIn = struct();
else
    unitsIn = sys.units;
end

lengthUnit = local_get_unit(unitsIn, 'length', 'angstrom');
alphaUnit  = local_get_unit(unitsIn, 'alpha',  'angstrom^3');
chargeUnit = local_get_unit(unitsIn, 'charge', 'elementary_charge');

lengthUnit = local_normalize_unit_label(lengthUnit);
alphaUnit  = local_normalize_unit_label(alphaUnit);
chargeUnit = local_normalize_unit_label(chargeUnit);

% -------------------------------------------------------------------------
% Positions
% -------------------------------------------------------------------------

if isfield(sys, 'site_pos') && ~isempty(sys.site_pos)
    if ~isnumeric(sys.site_pos) || size(sys.site_pos, 2) ~= 3
        error('io:convert_to_atomic_units:BadSitePos', ...
            'sys.site_pos must be a numeric N x 3 array.');
    end

    switch lengthUnit
        case 'angstrom'
            sys.site_pos = sys.site_pos * ANG2BOHR;

        case 'bohr'
            % already internal units

        otherwise
            error('io:convert_to_atomic_units:UnsupportedLengthUnit', ...
                'Unsupported length unit: %s', lengthUnit);
    end
end

% -------------------------------------------------------------------------
% Polarizabilities
% -------------------------------------------------------------------------

if isfield(sys, 'site_alpha') && ~isempty(sys.site_alpha)
    if ~isnumeric(sys.site_alpha)
        error('io:convert_to_atomic_units:BadSiteAlpha', ...
            'sys.site_alpha must be numeric.');
    end

    switch alphaUnit
        case 'angstrom^3'
            sys.site_alpha = sys.site_alpha / A3_PER_AU_ALPHA;

        case 'atomic_unit'
            % already internal units

        otherwise
            error('io:convert_to_atomic_units:UnsupportedAlphaUnit', ...
                'Unsupported polarizability unit: %s', alphaUnit);
    end
end

% -------------------------------------------------------------------------
% Charges
% -------------------------------------------------------------------------
% Charges stay in elementary-charge units. We only validate the label.

if isfield(sys, 'site_charge') && ~isempty(sys.site_charge)
    if ~isnumeric(sys.site_charge)
        error('io:convert_to_atomic_units:BadSiteCharge', ...
            'sys.site_charge must be numeric.');
    end
end

switch chargeUnit
    case 'elementary_charge'
        % already internal units

    otherwise
        error('io:convert_to_atomic_units:UnsupportedChargeUnit', ...
            'Unsupported charge unit: %s', chargeUnit);
end

% -------------------------------------------------------------------------
% Stamp canonical output unit metadata
% -------------------------------------------------------------------------

sys.units = struct();
sys.units.length = 'bohr';
sys.units.alpha = 'atomic_unit';
sys.units.charge = 'elementary_charge';

end

% =========================================================================
% Local helpers
% =========================================================================

function value = local_get_unit(unitsIn, fieldName, defaultValue)
if isfield(unitsIn, fieldName) && ~isempty(unitsIn.(fieldName))
    value = unitsIn.(fieldName);
else
    value = defaultValue;
end
end

function label = local_normalize_unit_label(label)
%LOCAL_NORMALIZE_UNIT_LABEL Normalize common unit label spellings.

if isstring(label)
    if numel(label) ~= 1
        error('io:convert_to_atomic_units:BadUnitLabel', ...
            'Unit labels must be scalar strings or character vectors.');
    end
    label = char(label);
end

if ~ischar(label)
    error('io:convert_to_atomic_units:BadUnitLabel', ...
        'Unit labels must be scalar strings or character vectors.');
end

label = lower(strtrim(label));

switch label
    case {'a', 'ang', 'angstrom', 'angstroms'}
        label = 'angstrom';

    case {'bohr', 'bohrs', 'au_length', 'a0'}
        label = 'bohr';

    case {'angstrom^3', 'angstrom3', 'ang^3', 'ang3', 'a^3', 'aa^3'}
        label = 'angstrom^3';

    case {'atomic_unit', 'atomic_units', 'au', 'a.u.', 'au_alpha', ...
            'bohr^3', 'bohr3'}
        label = 'atomic_unit';

    case {'elementary_charge', 'electron_charge', 'e', 'charge_e'}
        label = 'elementary_charge';

    otherwise
        % Keep the normalized spelling so the caller gets a useful error.
end
end