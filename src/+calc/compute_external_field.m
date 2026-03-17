function Eext = compute_external_field(sys, params)
%COMPUTE_EXTERNAL_FIELD External field at working sites from assigned charges.

fieldParams = struct();

if isfield(params, 'field') && ~isempty(params.field)
    fieldParams = params.field;
end

Eext = thole.induced_field_from_charges(sys, fieldParams);

end