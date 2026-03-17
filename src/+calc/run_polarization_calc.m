function result = run_polarization_calc(crystal, model, opts, params)
%RUN_POLARIZATION_CALC Main driver for one polarization calculation.
%
% Current implemented steps:
%   1) build working system
%   2) assign point-charge pattern if provided
%   3) compute external field from charges
%   4) build nonperiodic dipole interaction operator
%   5) solve induced dipoles
%   6) compute simple polarization energy
%
% Future steps:
%   - replace nonperiodic operator with periodic Ewald operator

if nargin < 4
    error('run_polarization_calc requires crystal, model, opts, and params.');
end

sys = builder.make_crystal_system(crystal, model, opts);

if isfield(model, 'charge_pattern') && ~isempty(model.charge_pattern)
    if isfield(sys, 'active_molecules') && ~isempty(sys.active_molecules)
        sys = builder.assign_point_charges(sys, model.charge_pattern);
    end
end

Eext = [];
if isfield(params, 'field') && isfield(params.field, 'include_external_charges')
    if params.field.include_external_charges
        Eext = calc.compute_external_field(sys, params);
    end
end

T = [];
mu = [];
scf = struct();
direct = struct();
Epol = [];

if ~isempty(Eext)
    T = ewald.assemble_nonperiodic_interaction_matrix(sys, params.ewald, params.scf);

    switch params.scf.solver
        case 'matrix_iterative'
            [mu, scf] = thole.solve_scf_matrix_iterative(sys, Eext, T, params.scf);

        case 'iterative'
            [mu, scf] = thole.solve_scf_iterative(sys, Eext, params.scf);

        case 'direct'
            [mu, direct] = thole.solve_scf_direct(sys, Eext, T);
            scf = struct();
            scf.converged = true;
            scf.nIter = 1;
            scf.history = 0;
            scf.used_matrix_solver = true;

        otherwise
            error('Unknown SCF solver: %s', params.scf.solver);
    end

    Epol = thole.dipole_energy(sys, mu, Eext);
end

result = struct();
result.sys = sys;
result.params = params;
result.Eext = Eext;
result.T = T;
result.mu = mu;
result.scf = scf;
result.direct = direct;
result.energy = struct();
result.energy.polarization = Epol;

result.status = 'direct-solver stage complete';
result.message = 'System built, charges assigned, operator assembled, induced dipoles solved.';

end