function params = default_params()
%DEFAULT_PARAMS Default numerical settings for polarization calculations.

params = struct();

params.ewald = struct();
params.ewald.mode = 'nonperiodic';   % 'nonperiodic' or 'periodic_triclinic'
params.ewald.alpha = 0.25;
params.ewald.rCut = 12.0;
params.ewald.kCut = 6;
params.ewald.boundary = 'tinfoil';
params.ewald.auto = false;
params.ewald.tol = 1e-8;
params.ewald.rcut_fraction = 0.9;

% Periodic real-space Thole correction controls
params.ewald.use_thole_real_space = false;
params.ewald.thole_a = [];

params.scf = struct();
params.scf.solver = 'gs';   % 'gs/sor', 'matrix_iterative', 'direct'
params.scf.tol = 1e-8;
params.scf.maxIter = 500;
params.scf.mixing = 0.5;
params.scf.initial_mu = [];
params.scf.softening = 0.0;
params.scf.use_thole = true;

params.field = struct();
params.field.include_external_charges = true;
params.field.exclude_self = true;
params.field.softening = 0.0;
params.field.target_mask = [];
params.field.source_mask = [];

params.output = struct();
params.output.verbose = true;
params.output.return_matrices = false;
params.output.return_fields = true;

end