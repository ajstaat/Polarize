function op = build_nonperiodic_dense_operator(sys, problem, opt)
%BUILD_NONPERIODIC_DENSE_OPERATOR Build dense nonperiodic polarization op.
%
% Private backend builder used by thole.make_polarization_operator.
%
% Returns:
%   op.kind    = 'dense_matrix'
%   op.backend = 'nonperiodic_paircache_dense'   for finite Rcut
%              = 'nonperiodic_allpairs_dense'    for Rcut = Inf

scfParams = struct();
scfParams.use_thole = opt.UseThole;
scfParams.softening = opt.Softening;
scfParams.rcut = opt.Rcut;
scfParams.use_mex = opt.UseMex;
scfParams.profile = opt.Profile;
scfParams.verbose = opt.Verbose;

[Tpol, info] = thole.assemble_nonperiodic_interaction_matrix(sys, problem, scfParams);

op = struct();
op.mode = 'nonperiodic';
op.kind = 'dense_matrix';

if isfinite(opt.Rcut)
    op.backend = 'nonperiodic_paircache_dense';
else
    op.backend = 'nonperiodic_allpairs_dense';
end

op.nPolSites = problem.nPolSites;
op.size = size(Tpol);

op.Tpol = Tpol;
op.apply = @(muVec) Tpol * muVec;

op.info = info;

op.capabilities = struct();
op.capabilities.apply = true;
op.capabilities.dense_matrix = true;
op.capabilities.row_update = false;

op.params = struct();
op.params.use_thole = opt.UseThole;
op.params.softening = opt.Softening;
op.params.rcut = opt.Rcut;
op.params.use_mex = opt.UseMex;

end