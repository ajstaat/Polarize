function Eext_vec = compute_external_field_vector(sys, params)
%COMPUTE_EXTERNAL_FIELD_VECTOR Compute stacked 3N x 1 external field vector.

Eext = calc.compute_external_field(sys, params);
Eext_vec = util.stack_xyz(Eext);

end