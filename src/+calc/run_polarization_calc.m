function result = run_polarization_calc(crystal, model, opts, params)
%RUN_POLARIZATION_CALC Main driver for one polarization calculation.
%
% Stages:
%   1) build crystal/supercell system
%   2) select active molecules / assign point charges
%   3) build external field from assigned charges
%   4) assemble dipole interaction operator
%   5) solve induced dipoles
%   6) optionally compute energy breakdown
%
% Supported SCF solver names:
%   'matrix_iterative', 'matrix', 'jacobi'
%   'iterative', 'fixed_point'
%   'gs', 'sor', 'gauss_seidel'
%   'direct'

    if nargin < 4
        error('run_polarization_calc requires crystal, model, opts, and params.');
    end

    % ---------------------------
    % Verbosity
    % ---------------------------
    verbose = false;
    if isfield(opts, 'verbose') && ~isempty(opts.verbose)
        verbose = logical(opts.verbose);
    elseif isfield(params, 'output') && isfield(params.output, 'verbose') ...
            && ~isempty(params.output.verbose)
        verbose = logical(params.output.verbose);
    end

        % 6) Requested outputs
        if computeEnergy || computeEnergyByMolecule || computeDipoleDipoleDecomp
            tEnergy = tic;
            if verbose
                fprintf('Computing requested post-SCF outputs...\n');
            end

            if computeEnergy
                energy = calc.compute_total_energy(sys, mu, Eext, T);
            else
                energy = struct();
            end

            if computeEnergyByMolecule
                energy_by_molecule = calc.compute_total_energy_by_molecule(sys, mu, Eext);
            else
                energy_by_molecule = table();
            end

            if computeDipoleDipoleDecomp
                dipole_dipole_decomp = calc.compute_dipole_dipole_decomposition(sys, mu, T);
            else
                dipole_dipole_decomp = struct();
            end

            if verbose
                fprintf('Post-SCF analysis finished in %.2f s\n', toc(tEnergy));
                if isfield(energy, 'total')
                    fprintf('Total polarization energy: %.10f\n', energy.total);
                end
            end
        else
            if verbose
                fprintf('Skipping post-SCF energy/decomposition outputs.\n');
            end
            energy = struct();
            energy_by_molecule = table();
            dipole_dipole_decomp = struct();
        end    

    % ---------------------------
    % 1) Build system
    % ---------------------------
    tBuild = tic;
    if verbose
        fprintf('Building crystal system...\n');
    end
    sys = builder.make_crystal_system(crystal, model, opts);
    if verbose
        nUniqueMol = numel(unique(sys.site_mol_id));
        fprintf('Built system in %.2f s: %d sites, %d unique molecule images\n', ...
            toc(tBuild), sys.n_sites, nUniqueMol);
    end

    % ---------------------------
    % 2) Assign charges
    % ---------------------------
    if isfield(sys, 'active_molecules') && ~isempty(sys.active_molecules)
        if verbose
            fprintf('Active molecules: %s\n', mat2str(sys.active_molecules));
        end

        tCharge = tic;
        if isfield(model, 'charge_patterns') && ~isempty(model.charge_patterns)
            if verbose
                fprintf('Assigning %d molecule-specific charge patterns...\n', ...
                    numel(model.charge_patterns));
            end
            sys = builder.assign_point_charges_multi(sys, model.charge_patterns);

        elseif isfield(model, 'charge_pattern') && ~isempty(model.charge_pattern)
            if verbose
                fprintf('Assigning shared charge pattern to active molecules...\n');
            end
            sys = builder.assign_point_charges(sys, model.charge_pattern);

        else
            if verbose
                fprintf('No charge pattern supplied.\n');
            end
            if ~isfield(sys, 'site_charge') || isempty(sys.site_charge)
                sys.site_charge = zeros(sys.n_sites, 1);
            end
        end

        if verbose
            fprintf('Charge assignment finished in %.2f s\n', toc(tCharge));
        end
    else
        if verbose
            fprintf('No active molecules selected.\n');
        end
        if ~isfield(sys, 'site_charge') || isempty(sys.site_charge)
            sys.site_charge = zeros(sys.n_sites, 1);
        end
    end

    % Optional: treat active charged molecules as fixed field sources only
    if isfield(opts, 'depolarizeActiveMolecules') && opts.depolarizeActiveMolecules
        if isfield(sys, 'active_molecules') && ~isempty(sys.active_molecules)
            if verbose
                fprintf('Depolarizing active molecules: %s\n', mat2str(sys.active_molecules));
            end
            sys = builder.depolarize_molecules(sys, sys.active_molecules);
        end
    end

    if verbose && isfield(sys, 'site_charge') && ~isempty(sys.site_charge)
        fprintf('Total assigned charge: %+0.6f e\n', sum(sys.site_charge));
        fprintf('Nonzero charged sites: %d\n', nnz(abs(sys.site_charge) > 1e-14));
        fprintf('Max |site charge|: %.6f e\n', max(abs(sys.site_charge)));
    end

    % ---------------------------
    % 3) Build external field only if there are nonzero charges
    % ---------------------------
    Eext = [];
    hasCharges = isfield(sys, 'site_charge') && ~isempty(sys.site_charge) ...
        && any(abs(sys.site_charge) > 1e-14);

    if hasCharges
        if isfield(params, 'field') && isfield(params.field, 'include_external_charges') ...
                && params.field.include_external_charges
            tField = tic;
            if verbose
                fprintf('Computing external field from assigned charges...\n');
            end
            Eext = calc.compute_external_field(sys, params);
            if verbose
                fprintf('External field computed in %.2f s\n', toc(tField));
            end
        else
            if verbose
                fprintf('External field disabled.\n');
            end
        end
    else
        if verbose
            fprintf('No nonzero assigned charges; skipping external-field build.\n');
        end
    end

    % Initialize outputs so result is always complete
    T = [];
    Tparts = struct();
    Tmeta = struct();
    mu = [];
    scf = struct();
    direct = struct();
    energy = struct();
    energy_by_molecule = table();
    dipole_dipole_decomp = struct();

    % ---------------------------
    % 4–6) Assemble operator, solve, compute requested outputs
    % ---------------------------
    if ~isempty(Eext)
        mode = 'nonperiodic';
        if isfield(params, 'ewald') && isfield(params.ewald, 'mode') && ~isempty(params.ewald.mode)
            mode = params.ewald.mode;
        end

        if verbose
            fprintf('Operator mode: %s\n', mode);
        end

        % 4) Assemble operator
        tOp = tic;
        switch lower(string(mode))
            case "nonperiodic"
                if verbose
                    fprintf('Assembling nonperiodic interaction matrix...\n');
                end
                T = ewald.assemble_nonperiodic_interaction_matrix(sys, params.ewald, params.scf);

            case "periodic_triclinic"
                if verbose
                    fprintf('Assembling periodic triclinic interaction matrix...\n');
                end
                ewaldParams = params.ewald;

                if isfield(ewaldParams, 'auto') && ewaldParams.auto
                    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
                        H = sys.super_lattice.';
                    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
                        H = sys.lattice.';
                    else
                        error('Need sys.super_lattice or sys.lattice for automatic Ewald parameter selection.');
                    end

                    if verbose
                        fprintf('Choosing automatic Ewald parameters...\n');
                    end

                    pauto = ewald.choose_ewald_params( ...
                        H, ewaldParams.tol, ewaldParams.boundary, ewaldParams.rcut_fraction);

                    ewaldParams.alpha = pauto.alpha;
                    ewaldParams.rcut = pauto.rcut;
                    ewaldParams.kcut = pauto.kcut;
                    ewaldParams.boundary = pauto.boundary;

                    if verbose
                        fprintf('Ewald params: alpha=%.4g, rcut=%.4g, kcut=%g\n', ...
                            ewaldParams.alpha, ewaldParams.rcut, ewaldParams.kcut);
                    end
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

        if verbose
            fprintf('Operator assembled in %.2f s\n', toc(tOp));
        end

        % 5) Solve induced dipoles
        tSolve = tic;
        if verbose
            fprintf('Starting SCF solve using solver "%s"...\n', params.scf.solver);
        end

        solver = lower(string(params.scf.solver));
        switch solver
            case {"matrix_iterative", "matrix", "jacobi"}
                [mu, scf] = thole.solve_scf_matrix_iterative(sys, Eext, T, params.scf);

            case {"iterative", "fixed_point"}
                [mu, scf] = thole.solve_scf_iterative(sys, Eext, params.scf);

            case {"gs", "sor", "gauss_seidel"}
                [mu, scf] = thole.solve_scf_sor(sys, Eext, T, params.scf);

            case {"direct"}
                [mu, direct] = thole.solve_scf_direct(sys, Eext, T);
                scf = struct();
                scf.method = 'direct';
                scf.converged = true;
                scf.nIter = 1;
                scf.history = 0;
                scf.used_matrix_solver = true;

            otherwise
                error('Unknown SCF solver: %s', params.scf.solver);
        end

        if verbose
            fprintf('SCF solve finished in %.2f s\n', toc(tSolve));
            if isfield(scf, 'method') && strcmpi(scf.method, 'direct')
                fprintf('Direct solve completed.\n');
            elseif isfield(scf, 'converged')
                fprintf('SCF converged: %d | iterations: %d\n', scf.converged, scf.nIter);
            end
        end

        % 6) Requested outputs
        if computeEnergy || computeEnergyByMolecule || computeDipoleDipoleDecomp
            tEnergy = tic;
            if verbose
                fprintf('Computing requested post-SCF outputs...\n');
            end

            if computeEnergy
                energy = calc.compute_total_energy(sys, mu, Eext, T);
            else
                energy = struct();
            end

            if computeEnergyByMolecule
                energy_by_molecule = calc.compute_total_energy_by_molecule(sys, mu, Eext);
            else
                energy_by_molecule = table();
            end

            if computeDipoleDipoleDecomp
                dipole_dipole_decomp = calc.compute_dipole_dipole_decomposition(sys, mu, T);
            else
                dipole_dipole_decomp = struct();
            end

            if verbose
                fprintf('Post-SCF analysis finished in %.2f s\n', toc(tEnergy));
                if isfield(energy, 'total')
                    fprintf('Total polarization energy: %.10f\n', energy.total);
                end
            end
        else
            if verbose
                fprintf('Skipping post-SCF energy/decomposition outputs.\n');
            end
            energy = struct();
            energy_by_molecule = table();
            dipole_dipole_decomp = struct();
        end

    else
        if verbose
            fprintf('Skipping operator assembly and SCF because Eext is empty.\n');
        end

        % Optional neutral defaults
        energy = struct();
        energy.polarization_self = 0;
        energy.external_charge_dipole = 0;
        energy.dipole_dipole = 0;
        energy.total = 0;
    end

    % ---------------------------
    % Package result
    % ---------------------------
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
    result.message = 'System built, operator assembled, induced dipoles solved, requested outputs reported.';

    if verbose
        fprintf('=== done ===\n\n');
    end
end