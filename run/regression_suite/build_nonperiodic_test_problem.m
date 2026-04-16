function [sys, polsys, Eext, problem, Tpol] = build_nonperiodic_test_problem(crystal, model, opts, params)
%BUILD_NONPERIODIC_TEST_PROBLEM Common regression-test setup.

    sys = builder.make_crystal_system(crystal, model, opts);
    sys = builder.assign_point_charges(sys, model.charge_pattern);

    if isfield(opts, 'depolarizeActiveMolecules') && opts.depolarizeActiveMolecules
        if isfield(sys, 'active_molecules') && ~isempty(sys.active_molecules)
            sys = builder.depolarize_molecules(sys, sys.active_molecules);
        end
    end

    polsys = builder.extract_polarization_system(sys, params);

    Eext = calc.compute_external_field(polsys, params);
    problem = thole.prepare_scf_problem(polsys, Eext, params.scf);
    [Tpol, ~] = ewald.assemble_nonperiodic_interaction_matrix(polsys, problem, params.ewald, params.scf);
end