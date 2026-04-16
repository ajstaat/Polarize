function [crystal, model, opts0, params0, molA, molB] = make_standard_nonperiodic_case(supercell_size)
%MAKE_STANDARD_NONPERIODIC_CASE Common tiny nonperiodic regression case.

    if nargin < 1 || isempty(supercell_size)
        supercell_size = [2 2 1];
    end

    crystal = struct();
    crystal.cellpar = [12.3, 14.1, 9.8, 90.0, 101.2, 88.4];
    crystal.lattice = [];
    crystal.frac_coords = [
        0.10 0.20 0.30
        0.15 0.22 0.35
        0.60 0.70 0.80
    ];
    crystal.cart_coords = [];
    crystal.mol_id = [1; 1; 2];
    crystal.site_label = {'C1'; 'H1'; 'N1'};
    crystal.site_type = {'C'; 'H'; 'N'};

    model = struct();
    model.polarizable_types = {'C', 'N'};
    model.alpha_by_type = struct();
    model.alpha_by_type.C = 10.0;
    model.alpha_by_type.N = 8.0;
    model.thole_a = 0.39;

    model.charge_pattern = struct();
    model.charge_pattern.site_label = {'C1', 'H1'};
    model.charge_pattern.delta_q = [0.10, -0.10];

    opts0 = struct();
    opts0.supercell_size = supercell_size;
    opts0.verbose = false;
    opts0.removeMolIDs = [];
    opts0.activeMolIDs = [];

    params0 = util.default_params();
    params0.scf.solver = 'matrix_iterative';
    params0.scf.tol = 1e-10;
    params0.scf.maxIter = 500;
    params0.scf.mixing = 0.35;
    params0.scf.use_thole = true;

    result0 = calc.run_polarization_calc(crystal, model, opts0, params0);
    sys0 = result0.sys;

    molA = builder.find_unique_molecule_id(sys0, 1, [0 0 0]);
    molB = builder.find_unique_molecule_id(sys0, 1, [1 0 0]);
end