function result = run_polarization_calc(crystal, model, opts, params)
%RUN_POLARIZATION_CALC Main driver for one polarization calculation.
%
% Stages:
%   1) build crystal/supercell system
%   2) select active molecules / assign charges
%   3) build external field from assigned charges
%   4) assemble dipole interaction operator (when required)
%   5) solve induced dipoles
%   6) optionally compute energy breakdown
%
% Supported SCF solver names:
%   Explicit-matrix solvers:
%       'matrix_iterative', 'matrix', 'jacobi'
%       'gs', 'sor', 'gauss_seidel'
%       'direct'
%
%   Matrix-free solver:
%       'iterative', 'fixed_point'
%
% Notes on current refactor state:
%   - nonperiodic mode uses active-space (polarizable-only) operator assembly
%     and the new problem-preparation interface for explicit-matrix solvers
%   - periodic_triclinic mode remains on the legacy full-space path for now
%   - builder.extract_polarization_system is the canonical boundary where
%     builder/raw units are converted into internal atomic units
%   - 'iterative' / 'fixed_point' is the matrix-free solver path and does
%     not assemble Tpol in nonperiodic mode

    if nargin < 4
        error('run_polarization_calc requires crystal, model, opts, and params.');
    end

    if nargin < 3 || isempty(opts)
        opts = struct();
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

    % ---------------------------
    % Requested outputs
    % ---------------------------
    computeEnergy = false;
    computeEnergyByMolecule = false;
    computeDipoleDipoleDecomp = false;

    if isfield(params, 'output') && ~isempty(params.output)
        if isfield(params.output, 'computeEnergy') && ~isempty(params.output.computeEnergy)
            computeEnergy = logical(params.output.computeEnergy);
        end
        if isfield(params.output, 'computeEnergyByMolecule') && ~isempty(params.output.computeEnergyByMolecule)
            computeEnergyByMolecule = logical(params.output.computeEnergyByMolecule);
        end
        if isfield(params.output, 'computeDipoleDipoleDecomp') && ~isempty(params.output.computeDipoleDipoleDecomp)
            computeDipoleDipoleDecomp = logical(params.output.computeDipoleDipoleDecomp);
        end
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
    % 2) Assign charges using refactored builder API
    % ---------------------------
    disablePolOnCharged = true;
    if isfield(opts, 'depolarizeActiveMolecules') && ~isempty(opts.depolarizeActiveMolecules)
        disablePolOnCharged = logical(opts.depolarizeActiveMolecules);
    end

    if isfield(sys, 'active_molecules') && ~isempty(sys.active_molecules)
        if verbose
            fprintf('Active molecules: %s\n', mat2str(sys.active_molecules));
        end

        tCharge = tic;

        if isfield(model, 'charge_files') && ~isempty(model.charge_files)
            if verbose
                fprintf('Assigning file-based charges to active molecules...\n');
            end

            requiredFileFields = {'template_coords', 'template_species'};
            for k = 1:numel(requiredFileFields)
                f = requiredFileFields{k};
                if ~isfield(model, f) || isempty(model.(f))
                    error('run_polarization_calc:MissingModelField', ...
                        'File-based charge assignment requires model.%s.', f);
                end
            end

            sys = builder.apply_molecule_charges( ...
                sys, sys.active_molecules, ...
                'Mode', 'file', ...
                'ChargeFiles', model.charge_files, ...
                'TemplateCoords', model.template_coords, ...
                'TemplateSpecies', model.template_species, ...
                'SetActive', true, ...
                'DisablePolarizabilityOnCharged', disablePolOnCharged, ...
                'ZeroExistingCharges', true, ...
                'Verbose', verbose);

        elseif isfield(model, 'total_charges') && ~isempty(model.total_charges)
            if verbose
                fprintf('Assigning uniform total charges to active molecules...\n');
            end

            sys = builder.apply_molecule_charges( ...
                sys, sys.active_molecules, ...
                'Mode', 'uniform', ...
                'TotalCharges', model.total_charges, ...
                'SetActive', true, ...
                'DisablePolarizabilityOnCharged', disablePolOnCharged, ...
                'ZeroExistingCharges', true, ...
                'Verbose', verbose);

        elseif isfield(model, 'charge_pattern') && ~isempty(model.charge_pattern)
            if isscalar(model.charge_pattern)
                if verbose
                    fprintf(['Assigning repeated uniform total charge %+g to ' ...
                             'all active molecules...\n'], model.charge_pattern);
                end

                totalCharges = repmat(model.charge_pattern, numel(sys.active_molecules), 1);

                sys = builder.apply_molecule_charges( ...
                    sys, sys.active_molecules, ...
                    'Mode', 'uniform', ...
                    'TotalCharges', totalCharges, ...
                    'SetActive', true, ...
                    'DisablePolarizabilityOnCharged', disablePolOnCharged, ...
                    'ZeroExistingCharges', true, ...
                    'Verbose', verbose);
            else
                error('run_polarization_calc:LegacyChargePatternUnsupported', ...
                    ['model.charge_pattern is non-scalar. In the refactored workflow, ' ...
                     'use model.total_charges for uniform mode or ' ...
                     'model.charge_files/model.template_coords/model.template_species ' ...
                     'for file mode.']);
            end

        elseif isfield(model, 'charge_patterns') && ~isempty(model.charge_patterns)
            error('run_polarization_calc:LegacyChargePatternsUnsupported', ...
                ['model.charge_patterns uses the older multi-pattern API. ' ...
                 'Please migrate this workflow to builder.apply_molecule_charges inputs.']);

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

    if verbose && isfield(sys, 'site_charge') && ~isempty(sys.site_charge)
        fprintf('Total assigned charge: %+0.6f e\n', sum(sys.site_charge));
        fprintf('Nonzero charged sites: %d\n', nnz(abs(sys.site_charge) > 1e-14));
        fprintf('Max |site charge|: %.6f e\n', max(abs(sys.site_charge)));
    end

    % ---------------------------
    % Canonical polarization system (converted to atomic units here)
    % ---------------------------
    polsys = builder.extract_polarization_system(sys, params);

    if verbose
        fprintf('Polarization system: %d sites (%d polarizable)\n', ...
            polsys.n_sites, nnz(polsys.site_is_polarizable));
        if isfield(polsys, 'units') && isstruct(polsys.units)
            fprintf('Internal units: length=%s, alpha=%s, charge=%s\n', ...
                polsys.units.length, polsys.units.alpha, polsys.units.charge);
        end
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
            Eext = calc.compute_external_field(polsys, params);
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

    % ---------------------------
    % Initialize outputs so result is always complete
    % ---------------------------
    T = [];
    Tpol = [];
    Tparts = struct();
    Tmeta = struct();
    opinfo = struct();
    problem = [];
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
        modeStr = lower(char(string(mode)));

        solver = lower(string(params.scf.solver));
        isMatrixFree = any(strcmp(solver, ["iterative","fixed_point"]));

        if verbose
            fprintf('Operator mode: %s\n', modeStr);
            if isMatrixFree
                fprintf('SCF path: matrix-free iterative solver\n');
            else
                fprintf('SCF path: explicit-matrix solver\n');
            end
        end

        % Prepare normalized SCF problem for explicit-matrix solvers and
        % downstream reporting. For matrix-free iterative, we do not require
        % Tpol, but we may still construct problem later if needed for
        % post-SCF reporting.
        if ~isMatrixFree || computeEnergy || computeEnergyByMolecule || computeDipoleDipoleDecomp
            problem = thole.prepare_scf_problem(polsys, Eext, params.scf);
        end

        % 4) Assemble operator when required
        tOp = tic;
        switch modeStr
            case 'nonperiodic'
                if ~isMatrixFree
                    if verbose
                        fprintf('Assembling nonperiodic interaction matrix...\n');
                    end

                    [Tpol, opinfo] = ewald.assemble_nonperiodic_interaction_matrix( ...
                        polsys, problem, params.ewald, params.scf);

                    T = [];
                    Tmeta = opinfo;

                    if verbose
                        fprintf('Operator assembled in %.2f s\n', toc(tOp));
                    end
                else
                    if verbose
                        fprintf('Skipping explicit Tpol assembly for matrix-free iterative solver.\n');
                    end
                end

            case 'periodic_triclinic'
                if isMatrixFree
                    error('Matrix-free iterative solver is not yet implemented for periodic_triclinic mode.');
                end

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
                        error(['Need sys.super_lattice or sys.lattice for automatic ' ...
                               'Ewald parameter selection.']);
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

                [T, Tparts, Tmeta] = ewald.assemble_periodic_interaction_matrix( ...
                    polsys, ewaldParams, params.scf);

                if verbose
                    fprintf('Operator assembled in %.2f s\n', toc(tOp));
                end

            otherwise
                error('Unknown operator mode: %s', modeStr);
        end

        % 5) Solve induced dipoles
        tSolve = tic;
        if verbose
            fprintf('Starting SCF solve using solver "%s"...\n', params.scf.solver);
        end

        switch solver
            case {"matrix_iterative", "matrix", "jacobi"}
                if strcmp(modeStr, 'nonperiodic')
                    [mu, scf] = thole.solve_scf_matrix_iterative(problem, Tpol);
                else
                    [mu, scf] = thole.solve_scf_matrix_iterative(polsys, Eext, T, params.scf);
                end

            case {"iterative", "fixed_point"}
                % Matrix-free iterative solver: apply dipole-field operator
                % directly without assembling dense Tpol.
                [mu, scf] = thole.solve_scf_iterative(polsys, Eext, params.scf);

            case {"gs", "sor", "gauss_seidel"}
                if strcmp(modeStr, 'nonperiodic')
                    [mu, scf] = thole.solve_scf_sor(problem, Tpol);
                else
                    [mu, scf] = thole.solve_scf_sor(polsys, Eext, T, params.scf);
                end

            case {"direct"}
                if strcmp(modeStr, 'nonperiodic')
                    [mu, direct] = thole.solve_scf_direct(problem, Tpol);
                else
                    [mu, direct] = thole.solve_scf_direct(polsys, Eext, T);
                end
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

        % For matrix-free nonperiodic solver, build problem/Tpol only if the
        % requested outputs require explicit-matrix post-SCF quantities.
        if strcmp(modeStr, 'nonperiodic') && isMatrixFree ...
                && (computeEnergy || computeEnergyByMolecule || computeDipoleDipoleDecomp)

            if isempty(problem)
                problem = thole.prepare_scf_problem(polsys, Eext, params.scf);
            end

            if isempty(Tpol)
                if verbose
                    fprintf('Assembling Tpol after matrix-free solve for requested post-SCF outputs...\n');
                end
                [Tpol, opinfo] = ewald.assemble_nonperiodic_interaction_matrix( ...
                    polsys, problem, params.ewald, params.scf);
                Tmeta = opinfo;
            end
        end

        % 6) Requested outputs
        if computeEnergy || computeEnergyByMolecule || computeDipoleDipoleDecomp
            tEnergy = tic;
            if verbose
                fprintf('Computing requested post-SCF outputs...\n');
            end

            if computeEnergy
                if strcmp(modeStr, 'nonperiodic')
                    if isempty(problem) || isempty(Tpol)
                        error(['Nonperiodic energy reporting currently requires problem and Tpol. ' ...
                               'These were not available.']);
                    end
                    energy = calc.compute_total_energy_active_space(polsys, problem, mu, Eext, Tpol);
                else
                    energy = calc.compute_total_energy(polsys, mu, Eext, T);
                end
            else
                energy = struct();
            end

            if computeEnergyByMolecule
                if strcmp(modeStr, 'nonperiodic')
                    if isempty(problem) || isempty(Tpol)
                        error(['Nonperiodic per-molecule energy reporting currently requires problem and Tpol. ' ...
                               'These were not available.']);
                    end
                    energy_by_molecule = calc.compute_total_energy_by_molecule_active_space( ...
                        polsys, problem, mu, Eext, Tpol);
                else
                    energy_by_molecule = calc.compute_total_energy_by_molecule(polsys, mu, Eext);
                end
            else
                energy_by_molecule = table();
            end

            if computeDipoleDipoleDecomp
                if strcmp(modeStr, 'nonperiodic')
                    if isempty(problem) || isempty(Tpol)
                        error(['Nonperiodic dipole-dipole decomposition currently requires problem and Tpol. ' ...
                               'These were not available.']);
                    end
                    dipole_dipole_decomp = calc.compute_dipole_dipole_decomposition_active_space( ...
                        polsys, problem, mu, Tpol);
                else
                    dipole_dipole_decomp = calc.compute_dipole_dipole_decomposition(polsys, mu, T);
                end
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
    result.polsys = polsys;
    result.params = params;
    result.Eext = Eext;
    result.T = T;
    result.Tpol = Tpol;
    result.Tparts = Tparts;
    result.Tmeta = Tmeta;
    result.problem = problem;
    result.mu = mu;
    result.scf = scf;
    result.direct = direct;
    result.energy = energy;
    result.energy_by_molecule = energy_by_molecule;
    result.dipole_dipole_decomposition = dipole_dipole_decomp;
    result.status = 'energy-reporting stage complete';
    result.message = ['System built, operator assembled or applied, induced dipoles solved, ' ...
                      'requested outputs reported.'];

    if verbose
        fprintf('=== done ===\n\n');
    end
end