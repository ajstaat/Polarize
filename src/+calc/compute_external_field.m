function Eext = compute_external_field(sys, params)
%COMPUTE_EXTERNAL_FIELD External field at working sites from assigned charges.
%
% Input
%   sys     canonical polarization-system struct in atomic units
%   params  calculation params; params.field is optional
%
% Output
%   Eext    external electric field at working sites

    io.assert_atomic_units(sys);

    fieldParams = struct();
    if nargin >= 2 && isfield(params, 'field') && ~isempty(params.field)
        fieldParams = params.field;
    end

    Eext = thole.induced_field_from_charges(sys, fieldParams);
end