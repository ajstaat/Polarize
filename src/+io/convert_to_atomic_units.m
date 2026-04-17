function sys = convert_to_atomic_units(sys)
%CONVERT_TO_ATOMIC_UNITS Convert system geometry/polarizabilities to a.u.
%
% Internal convention after conversion:
%   - positions: bohr
%   - polarizabilities: atomic units
%   - charges: elementary charge (unchanged)
%
% Required input fields:
%   sys.site_pos    [N x 3]
%   sys.site_alpha  [N x 1]
%
% Optional metadata:
%   sys.units.length   e.g. 'angstrom', 'bohr'
%   sys.units.alpha    e.g. 'angstrom^3', 'a.u.'
%   sys.units.charge   e.g. 'e'
%
% Behavior:
%   - If units metadata is present, convert based on it.
%   - If units metadata is absent, assume legacy crystal-builder inputs:
%       positions in Angstrom
%       polarizabilities in Angstrom^3
%   - Charges are left unchanged.

    if ~isstruct(sys)
        error('sys must be a struct.');
    end

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('sys.site_pos is missing or empty.');
    end
    if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
        error('sys.site_alpha is missing or empty.');
    end

    if size(sys.site_pos, 2) ~= 3
        error('sys.site_pos must be N x 3.');
    end

    nSites = size(sys.site_pos, 1);
    if numel(sys.site_alpha) ~= nSites
        error('sys.site_alpha must have length N.');
    end

    ANG2BOHR  = 1.8897259886;
    A3_PER_AU = 0.148184711;

    % ---------------------------------------------------------------------
    % Determine current units
    % ---------------------------------------------------------------------
    lengthUnit = 'angstrom';
    alphaUnit  = 'angstrom^3';
    chargeUnit = 'elementary_charge';

    if isfield(sys, 'units') && ~isempty(sys.units)
        if isfield(sys.units, 'length') && ~isempty(sys.units.length)
            lengthUnit = normalize_length_unit(sys.units.length);
        end
        if isfield(sys.units, 'alpha') && ~isempty(sys.units.alpha)
            alphaUnit = normalize_alpha_unit(sys.units.alpha);
        end
        if isfield(sys.units, 'charge') && ~isempty(sys.units.charge)
            chargeUnit = normalize_charge_unit(sys.units.charge);
        end
    end

    % ---------------------------------------------------------------------
    % Convert positions
    % ---------------------------------------------------------------------
    switch lengthUnit
        case 'bohr'
            % already in internal units
        case 'angstrom'
            sys.site_pos = sys.site_pos * ANG2BOHR;
        otherwise
            error('Unsupported length unit: %s', lengthUnit);
    end

    % ---------------------------------------------------------------------
    % Convert polarizabilities
    % ---------------------------------------------------------------------
    sys.site_alpha = sys.site_alpha(:);
    switch alphaUnit
        case 'atomic_unit'
            % already in internal units
        case 'angstrom^3'
            sys.site_alpha = sys.site_alpha / A3_PER_AU;
        otherwise
            error('Unsupported polarizability unit: %s', alphaUnit);
    end

    % ---------------------------------------------------------------------
    % Charges
    % ---------------------------------------------------------------------
    % Charges stay in elementary-charge units. We only validate the label.
    switch chargeUnit
        case 'elementary_charge'
            % unchanged
        otherwise
            error('Unsupported charge unit: %s', chargeUnit);
    end

    % ---------------------------------------------------------------------
    % Stamp normalized metadata
    % ---------------------------------------------------------------------
    if ~isfield(sys, 'units') || isempty(sys.units)
        sys.units = struct();
    end

    sys.units.length = 'bohr';
    sys.units.alpha  = 'atomic_unit';
    sys.units.charge = 'elementary_charge';

    if ~isfield(sys, 'units_history') || isempty(sys.units_history)
        sys.units_history = {};
    end

    sys.units_history{end+1} = struct( ...
        'action', 'convert_to_atomic_units', ...
        'from_length', lengthUnit, ...
        'from_alpha', alphaUnit, ...
        'from_charge', chargeUnit, ...
        'to_length', 'bohr', ...
        'to_alpha', 'atomic_unit', ...
        'to_charge', 'elementary_charge');
end

% =========================================================================
% Unit normalization helpers
% =========================================================================

function u = normalize_length_unit(uin)
    s = lower(strtrim(string(uin)));
    if any(strcmp(s, ["angstrom","ang","angs","å","a"]))
        u = 'angstrom';
    elseif any(strcmp(s, ["bohr","au","a.u.","atomic_unit","atomic units"]))
        u = 'bohr';
    else
        error('Unrecognized length unit: %s', char(s));
    end
end

function u = normalize_alpha_unit(uin)
    s = lower(strtrim(string(uin)));
    if any(strcmp(s, ["angstrom^3","a^3","ang^3","angs^3"]))
        u = 'angstrom^3';
    elseif any(strcmp(s, ["atomic_unit","atomic units","au","a.u."]))
        u = 'atomic_unit';
    else
        error('Unrecognized polarizability unit: %s', char(s));
    end
end

function u = normalize_charge_unit(uin)
    s = lower(strtrim(string(uin)));
    if any(strcmp(s, ["e","electron_charge","elementary_charge"]))
        u = 'elementary_charge';
    else
        error('Unrecognized charge unit: %s', char(s));
    end
end