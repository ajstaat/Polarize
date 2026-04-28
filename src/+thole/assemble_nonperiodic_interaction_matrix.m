function [Tpol, opinfo] = assemble_nonperiodic_interaction_matrix(sys, problem, scfParams)
%ASSEMBLE_NONPERIODIC_INTERACTION_MATRIX Build dense nonperiodic Thole operator.
%
% [Tpol, opinfo] = thole.assemble_nonperiodic_interaction_matrix(sys, problem)
% [Tpol, opinfo] = thole.assemble_nonperiodic_interaction_matrix(sys, problem, scfParams)
%
% Builds the active-space dipole-dipole interaction matrix Tpol such that:
%
%   E_dip_pol_vec = Tpol * mu_pol_vec
%
% where both vectors are stacked over the polarizable active sites in the
% ordering defined by problem.activeSites.
%
% Inputs
%   sys
%       polarization system in atomic units:
%           sys.site_pos              N x 3, bohr
%           sys.site_alpha            N x 1, atomic units
%           sys.site_is_polarizable   N x 1 logical
%           sys.thole_a               scalar
%
%   problem
%       output from thole.prepare_scf_problem, with fields:
%           problem.activeSites
%           problem.nPolSites
%
%   scfParams optional struct fields:
%       .use_thole   logical, default true
%       .softening   scalar, default 0
%       .rcut        scalar cutoff in bohr, default Inf
%       .verbose     logical, default false
%
% Output
%   Tpol
%       3*nPolSites x 3*nPolSites dense active-space interaction matrix
%
%   opinfo
%       diagnostic struct:
%           .nPolSites
%           .nPairBlocks
%           .nPairBlocksKept
%           .nPairBlocksSkippedCutoff
%           .use_thole
%           .softening
%           .rcut
%           .use_cutoff
%           .assembly_time
%
% Notes
%   - Self blocks are zero.
%   - Pair blocks are filled symmetrically.
%   - This is the simple dense nonperiodic path. It does not use spatial
%     caches, MEX, or periodic/Ewald machinery.

if nargin < 3 || isempty(scfParams)
    scfParams = struct();
end

validate_sys(sys);
validate_problem(sys, problem);

useThole = local_get_field(scfParams, 'use_thole', true);
softening = local_get_field(scfParams, 'softening', 0.0);
rcut = local_get_field(scfParams, 'rcut', Inf);
verbose = local_get_field(scfParams, 'verbose', false);

if ~(islogical(useThole) && isscalar(useThole))
    error('thole:assemble_nonperiodic_interaction_matrix:BadUseThole', ...
        'scfParams.use_thole must be a logical scalar.');
end

if ~(isnumeric(softening) && isscalar(softening) && isfinite(softening) && softening >= 0)
    error('thole:assemble_nonperiodic_interaction_matrix:BadSoftening', ...
        'scfParams.softening must be a finite nonnegative scalar.');
end

if ~(isnumeric(rcut) && isscalar(rcut) && rcut > 0)
    error('thole:assemble_nonperiodic_interaction_matrix:BadRcut', ...
        'scfParams.rcut must be a positive scalar or Inf.');
end

if ~(islogical(verbose) && isscalar(verbose))
    error('thole:assemble_nonperiodic_interaction_matrix:BadVerbose', ...
        'scfParams.verbose must be a logical scalar.');
end

io.assert_atomic_units(sys);

sites = problem.activeSites(:);
nPol = problem.nPolSites;

Tpol = zeros(3*nPol, 3*nPol);

opts = struct();
opts.use_thole = useThole;
opts.softening = softening;

useCutoff = isfinite(rcut);
rcut2 = rcut^2;

nPairBlocks = nPol * (nPol - 1) / 2;
nPairBlocksKept = 0;
nPairBlocksSkippedCutoff = 0;

tStart = tic;

for a = 1:(nPol - 1)
    i = sites(a);
    ri = sys.site_pos(i, :);
    alpha_i = sys.site_alpha(i);

    Ia = local_block_indices(a);

    for b = (a + 1):nPol
        j = sites(b);
        rj = sys.site_pos(j, :);
        alpha_j = sys.site_alpha(j);

        dr = ri - rj;
        r2 = sum(dr.^2);

        if useCutoff && r2 > rcut2
            nPairBlocksSkippedCutoff = nPairBlocksSkippedCutoff + 1;
            continue;
        end

        Ib = local_block_indices(b);

        Tij = thole.dipole_tensor_block( ...
            ri, rj, alpha_i, alpha_j, sys.thole_a, opts);

        Tpol(Ia, Ib) = Tij;
        Tpol(Ib, Ia) = Tij.';

        nPairBlocksKept = nPairBlocksKept + 1;
    end
end

assemblyTime = toc(tStart);

opinfo = struct();
opinfo.nPolSites = nPol;
opinfo.nPairBlocks = nPairBlocks;
opinfo.nPairBlocksKept = nPairBlocksKept;
opinfo.nPairBlocksSkippedCutoff = nPairBlocksSkippedCutoff;
opinfo.use_thole = useThole;
opinfo.softening = softening;
opinfo.rcut = rcut;
opinfo.use_cutoff = useCutoff;
opinfo.assembly_time = assemblyTime;

if verbose
    fprintf('assemble_nonperiodic_interaction_matrix:\n');
    fprintf('  nPolSites                 = %d\n', nPol);
    fprintf('  Tpol size                 = %d x %d\n', size(Tpol,1), size(Tpol,2));
    fprintf('  pair blocks total         = %d\n', nPairBlocks);
    fprintf('  pair blocks kept          = %d\n', nPairBlocksKept);
    fprintf('  pair blocks skipped cutoff= %d\n', nPairBlocksSkippedCutoff);
    fprintf('  use_thole                 = %d\n', useThole);
    fprintf('  softening                 = %.6g\n', softening);
    fprintf('  rcut                      = %.6g\n', rcut);
    fprintf('  assembly time             = %.6f s\n', assemblyTime);
end

end

% =========================================================================
% Validation helpers
% =========================================================================

function validate_sys(sys)

required = {
    'site_pos'
    'site_alpha'
    'site_is_polarizable'
    'thole_a'
    'units'
};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(sys, name) || isempty(sys.(name))
        error('thole:assemble_nonperiodic_interaction_matrix:MissingField', ...
            'sys.%s is required and missing/empty.', name);
    end
end

if ~isnumeric(sys.site_pos) || size(sys.site_pos, 2) ~= 3
    error('thole:assemble_nonperiodic_interaction_matrix:BadSitePos', ...
        'sys.site_pos must be N x 3 numeric.');
end

nSites = size(sys.site_pos, 1);

if numel(sys.site_alpha) ~= nSites
    error('thole:assemble_nonperiodic_interaction_matrix:BadAlphaLength', ...
        'sys.site_alpha must have one entry per site.');
end

if numel(sys.site_is_polarizable) ~= nSites
    error('thole:assemble_nonperiodic_interaction_matrix:BadMaskLength', ...
        'sys.site_is_polarizable must have one entry per site.');
end

if ~(isnumeric(sys.thole_a) && isscalar(sys.thole_a) && isfinite(sys.thole_a))
    error('thole:assemble_nonperiodic_interaction_matrix:BadTholeA', ...
        'sys.thole_a must be a finite scalar.');
end

end

function validate_problem(sys, problem)

if ~isstruct(problem)
    error('thole:assemble_nonperiodic_interaction_matrix:BadProblem', ...
        'problem must be a struct from thole.prepare_scf_problem.');
end

required = {'activeSites', 'nPolSites'};

for k = 1:numel(required)
    name = required{k};

    if ~isfield(problem, name) || isempty(problem.(name))
        error('thole:assemble_nonperiodic_interaction_matrix:MissingProblemField', ...
            'problem.%s is required and missing/empty.', name);
    end
end

sites = problem.activeSites(:);
nPol = problem.nPolSites;

if ~(isnumeric(nPol) && isscalar(nPol) && nPol == numel(sites))
    error('thole:assemble_nonperiodic_interaction_matrix:BadNPolSites', ...
        'problem.nPolSites must equal numel(problem.activeSites).');
end

nSites = size(sys.site_pos, 1);

if any(sites < 1) || any(sites > nSites) || any(sites ~= round(sites))
    error('thole:assemble_nonperiodic_interaction_matrix:BadActiveSites', ...
        'problem.activeSites contains invalid site indices.');
end

if any(~logical(sys.site_is_polarizable(sites)))
    error('thole:assemble_nonperiodic_interaction_matrix:NonpolarizableActiveSite', ...
        'All problem.activeSites must be polarizable in sys.site_is_polarizable.');
end

end

% =========================================================================
% Small helpers
% =========================================================================

function idx = local_block_indices(k)

idx = (3*(k-1) + 1):(3*k);

end

function value = local_get_field(s, name, defaultValue)

if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end

end