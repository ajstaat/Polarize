function Eext = compute_external_field(sys, params)
%COMPUTE_EXTERNAL_FIELD External field at working sites from assigned charges.
%
% Input
%   sys     canonical polarization-system struct in atomic units
%   params  calculation params; params.field is optional
%
% Output
%   Eext    external electric field at working sites
%
% Behavior
%   fieldParams.mode:
%       'nonperiodic'  -> thole.induced_field_from_charges
%       'periodic'     -> thole.induced_field_from_charges_periodic
%
%   If params.field.use_thole_damping is not explicitly provided, it will
%   inherit from either:
%       params.use_thole
%   or  params.scf.use_thole
%   if available.

    io.assert_atomic_units(sys);

    fieldParams = struct();
    if nargin >= 2 && isfield(params, 'field') && ~isempty(params.field)
        fieldParams = params.field;
    end

    % Default mode
    if ~isfield(fieldParams, 'mode') || isempty(fieldParams.mode)
        fieldParams.mode = 'nonperiodic';
    end

    % Inherit top-level use_thole unless explicitly overridden for field.
    if ~isfield(fieldParams, 'use_thole_damping') || isempty(fieldParams.use_thole_damping)
        inheritedUseThole = [];

        if nargin >= 2 && isfield(params, 'use_thole') && ~isempty(params.use_thole)
            inheritedUseThole = logical(params.use_thole);
        elseif nargin >= 2 && isfield(params, 'scf') && ~isempty(params.scf) && ...
               isfield(params.scf, 'use_thole') && ~isempty(params.scf.use_thole)
            inheritedUseThole = logical(params.scf.use_thole);
        end

        if ~isempty(inheritedUseThole)
            fieldParams.use_thole_damping = inheritedUseThole;
        end
    end

    modeStr = lower(string(fieldParams.mode));

    switch modeStr
        case "nonperiodic"
            Eext = thole.induced_field_from_charges(sys, fieldParams);

        case "periodic"
            Eext = thole.induced_field_from_charges_periodic(sys, fieldParams);

        otherwise
            error('calc:compute_external_field:UnknownMode', ...
                'Unknown field mode "%s". Use "nonperiodic" or "periodic".', modeStr);
    end
end