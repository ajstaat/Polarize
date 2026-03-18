function result = run_polarization_calc(crystal, model, opts, params)
%RUN_POLARIZATION_CALC Main driver for one polarization calculation.

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
energy = struct();
energy_by_molecule = table();
dipole_dipole_decomp = struct();

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

    energy = calc.compute_total_energy(sys, mu, Eext, T);
    energy_by_molecule = calc.compute_total_energy_by_molecule(sys, mu, Eext);
    dipole_dipole_decomp = calc.compute_dipole_dipole_decomposition(sys, mu, T);
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
result.energy = energy;
result.energy_by_molecule = energy_by_molecule;
result.dipole_dipole_decomposition = dipole_dipole_decomp;

result.status = 'energy-reporting stage complete';
result.message = 'System built, operator assembled, induced dipoles solved, energy breakdown reported.';

end