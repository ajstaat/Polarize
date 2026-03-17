function result = run_polarization_calc(crystal, model, opts, params)
%RUN_POLARIZATION_CALC Main driver for one polarization calculation.
%
% Current implemented steps:
%   1) build working system
%   2) assign point-charge pattern if provided
%   3) compute external field from charges
%   4) build dipole interaction operator (nonperiodic or periodic)
%   5) solve induced dipoles
%   6) compute simple polarization energy
%
% Notes
%   - periodic_triclinic currently uses the Ewald dipole matrix operator
%   - nonperiodic path can still use Thole damping in the SCF kernel
%   - direct and matrix_iterative solvers use the explicit operator T

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
Tparts = struct();
Tmeta = struct();
mu = [];
scf = struct();
direct = struct();
Epol = [];

if ~isempty(Eext)
    mode = 'nonperiodic';
    if isfield(params, 'ewald') && isfield(params.ewald, 'mode') && ~isempty(params.ewald.mode)
        mode = params.ewald.mode;
    end

    switch mode
        case 'nonperiodic'
            T = ewald.assemble_nonperiodic_interaction_matrix(sys, params.ewald, params.scf);

        case 'periodic_triclinic'
            ewaldParams = params.ewald;

            % Optional automatic parameter selection from the system lattice
            if isfield(ewaldParams, 'auto') && ewaldParams.auto
                if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
                    H = sys.super_lattice.';
                elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
                    H = sys.lattice.';
                else
                    error('Need sys.super_lattice or sys.lattice for automatic Ewald parameter selection.');
                end

                pauto = ewald.choose_ewald_params( ...
                    H, ewaldParams.tol, ewaldParams.boundary, ewaldParams.rcut_fraction);

                ewaldParams.alpha = pauto.alpha;
                ewaldParams.rcut  = pauto.rcut;
                ewaldParams.kcut  = pauto.kcut;
                ewaldParams.boundary = pauto.boundary;
            else
                % allow older field names rCut/kCut as aliases
                if ~isfield(ewaldParams, 'rcut') && isfield(ewaldParams, 'rCut')
                    ewaldParams.rcut = ewaldParams.rCut;
                end
                if ~isfield(ewaldParams, 'kcut') && isfield(ewaldParams, 'kCut')
                    ewaldParams.kcut = ewaldParams.kCut;
                end
            end

            [T, Tparts, Tmeta] = ewald.assemble_periodic_interaction_matrix(sys, ewaldParams, params.scf);

        otherwise
            error('Unknown operator mode: %s', mode);
    end

    switch params.scf.solver
        case 'matrix_iterative'
            [mu, scf] = thole.solve_scf_matrix_iterative(sys, Eext, T, params.scf);

        case 'iterative'
            % This path uses the field-by-field kernel rather than explicit T.
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
result.Tparts = Tparts;
result.Tmeta = Tmeta;
result.mu = mu;
result.scf = scf;
result.direct = direct;
result.energy = struct();
result.energy.polarization = Epol;

result.status = 'operator-integrated stage complete';
result.message = 'System built, charges assigned, operator assembled, induced dipoles solved.';

end