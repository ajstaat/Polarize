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
% Finite cutoff path:
%   uses geom.build_nonperiodic_pair_cache with site_mask restricted to
%   problem.activeSites. That cache is full-site indexed, so this assembler
%   maps full site indices back to active-space block indices.
%
% Infinite cutoff path:
%   uses a simple all-pairs loop over problem.activeSites.
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
%       .use_mex     logical, default true
%       .profile     logical, default false
%       .verbose     logical, default false
%
% Output
%   Tpol
%       3*nPolSites x 3*nPolSites dense active-space interaction matrix
%
%   opinfo
%       diagnostic struct

if nargin < 3 || isempty(scfParams)
    scfParams = struct();
end

validate_sys(sys);
validate_problem(sys, problem);

useThole = local_get_field(scfParams, 'use_thole', true);
softening = local_get_field(scfParams, 'softening', 0.0);
rcut = local_get_field(scfParams, 'rcut', Inf);
useMex = local_get_field(scfParams, 'use_mex', true);
profile = local_get_field(scfParams, 'profile', false);
verbose = local_get_field(scfParams, 'verbose', false);

validate_options(useThole, softening, rcut, useMex, profile, verbose);

io.assert_atomic_units(sys);

sites = problem.activeSites(:);
nPol = problem.nPolSites;

Tpol = zeros(3*nPol, 3*nPol);

opts = struct();
opts.use_thole = useThole;
opts.softening = softening;

useCutoff = isfinite(rcut);
nPairBlocks = nPol * (nPol - 1) / 2;

tStart = tic;

if useCutoff
    [Tpol, fillInfo] = local_fill_from_pair_cache(Tpol, sys, sites, rcut, opts, ...
        useMex, profile, verbose);
    assemblyBackend = 'pair_cache';
else
    [Tpol, fillInfo] = local_fill_all_pairs(Tpol, sys, sites, opts);
    assemblyBackend = 'all_pairs';
end

assemblyTime = toc(tStart);

opinfo = struct();
opinfo.nPolSites = nPol;
opinfo.nPairBlocks = nPairBlocks;
opinfo.nPairBlocksKept = fillInfo.nPairBlocksKept;
opinfo.nPairBlocksSkippedCutoff = fillInfo.nPairBlocksSkippedCutoff;
opinfo.use_thole = useThole;
opinfo.softening = softening;
opinfo.rcut = rcut;
opinfo.use_cutoff = useCutoff;
opinfo.use_mex = useMex;
opinfo.profile = profile;
opinfo.assembly_backend = assemblyBackend;
opinfo.assembly_time = assemblyTime;

if isfield(fillInfo, 'cache_info')
    opinfo.cache_info = fillInfo.cache_info;
end

if isfield(fillInfo, 'cache_time')
    opinfo.cache_time = fillInfo.cache_time;
end

if isfield(fillInfo, 'fill_time')
    opinfo.fill_time = fillInfo.fill_time;
end

if verbose
    fprintf('assemble_nonperiodic_interaction_matrix:\n');
    fprintf('  nPolSites                 = %d\n', nPol);
    fprintf('  Tpol size                 = %d x %d\n', size(Tpol,1), size(Tpol,2));
    fprintf('  assembly backend          = %s\n', assemblyBackend);
    fprintf('  pair blocks total         = %d\n', nPairBlocks);
    fprintf('  pair blocks kept          = %d\n', fillInfo.nPairBlocksKept);
    fprintf('  pair blocks skipped cutoff= %d\n', fillInfo.nPairBlocksSkippedCutoff);
    fprintf('  use_thole                 = %d\n', useThole);
    fprintf('  softening                 = %.6g\n', softening);
    fprintf('  rcut                      = %.6g\n', rcut);
    fprintf('  use_mex                   = %d\n', useMex);
    fprintf('  assembly time             = %.6f s\n', assemblyTime);

    if isfield(opinfo, 'cache_time')
        fprintf('  cache time                = %.6f s\n', opinfo.cache_time);
    end

    if isfield(opinfo, 'fill_time')
        fprintf('  fill time                 = %.6f s\n', opinfo.fill_time);
    end
end

end

% =========================================================================
% Assembly paths
% =========================================================================

function [Tpol, info] = local_fill_from_pair_cache(Tpol, sys, activeSites, rcut, tensorOpts, ...
    useMex, profile, verbose)

nSites = size(sys.site_pos, 1);
nPol = numel(activeSites);

siteMask = false(nSites, 1);
siteMask(activeSites) = true;

cacheOpts = struct();
cacheOpts.rcut = rcut;
cacheOpts.site_mask = siteMask;
cacheOpts.use_mex = useMex;
cacheOpts.profile = profile || verbose;

tCache = tic;
cache = geom.build_nonperiodic_pair_cache(sys, cacheOpts);
cacheTime = toc(tCache);

fullToActive = zeros(nSites, 1);
fullToActive(activeSites) = 1:nPol;

pair_i = cache.pair_i(:);
pair_j = cache.pair_j(:);

a = fullToActive(pair_i);
b = fullToActive(pair_j);

if any(a == 0) || any(b == 0)
    error('thole:assemble_nonperiodic_interaction_matrix:BadPairCacheMask', ...
        'Pair cache returned a site outside the active-site mask.');
end

if any(a == b)
    keep = (a ~= b);
    pair_i = pair_i(keep);
    pair_j = pair_j(keep);
    a = a(keep);
    b = b(keep);

    cache.dr = cache.dr(keep, :);
    cache.r_bare = cache.r_bare(keep);
    cache.r2_bare = cache.r2_bare(keep);
    cache.inv_r3_bare = cache.inv_r3_bare(keep);
    cache.inv_r5_bare = cache.inv_r5_bare(keep);

    if isfield(cache, 'thole_f3')
        cache.thole_f3 = cache.thole_f3(keep);
        cache.thole_f5 = cache.thole_f5(keep);
    end
end

tFill = tic;

% -------------------------------------------------------------------------
% Vectorized tensor construction.
%
% dipole_tensor_block uses rvec = ri - rj.
% build_nonperiodic_pair_cache stores dr = rj - ri in the fallback path.
% The tensor depends on rvec*rvec', so the sign does not matter.
% -------------------------------------------------------------------------

dx = cache.dr(:, 1);
dy = cache.dr(:, 2);
dz = cache.dr(:, 3);

if tensorOpts.softening == 0
    invR3 = cache.inv_r3_bare(:);
    invR5 = cache.inv_r5_bare(:);
else
    r2 = cache.r2_bare(:) + tensorOpts.softening^2;
    r = sqrt(r2);
    invR3 = 1 ./ (r2 .* r);
    invR5 = invR3 ./ r2;
end

if tensorOpts.use_thole
    if isfield(cache, 'thole_f3') && isfield(cache, 'thole_f5')
        f3 = cache.thole_f3(:);
        f5 = cache.thole_f5(:);
    else
        tf = thole.thole_f3f5_factors( ...
            cache.r_bare(:), ...
            sys.site_alpha(pair_i), ...
            sys.site_alpha(pair_j), ...
            sys.thole_a);
        f3 = tf.f3(:);
        f5 = tf.f5(:);
    end
else
    f3 = ones(size(invR3));
    f5 = ones(size(invR5));
end

c5 = 3 .* f5 .* invR5;
c3 = f3 .* invR3;

% Tensor components for pair (i,j):
%
% T = c5 * dr*dr' - c3 * I
%
% Fill T(a,b) and T(b,a). The tensor is symmetric, so transpose has the
% same scalar components, but we still assign both active-space blocks.
Txx = c5 .* dx .* dx - c3;
Txy = c5 .* dx .* dy;
Txz = c5 .* dx .* dz;

Tyy = c5 .* dy .* dy - c3;
Tyz = c5 .* dy .* dz;

Tzz = c5 .* dz .* dz - c3;

n = size(Tpol, 1);

ra1 = 3*(a - 1) + 1;
ra2 = ra1 + 1;
ra3 = ra1 + 2;

rb1 = 3*(b - 1) + 1;
rb2 = rb1 + 1;
rb3 = rb1 + 2;

% Upper/off-diagonal block: rows a, cols b
idx = [
    sub2ind([n n], ra1, rb1)
    sub2ind([n n], ra1, rb2)
    sub2ind([n n], ra1, rb3)
    sub2ind([n n], ra2, rb1)
    sub2ind([n n], ra2, rb2)
    sub2ind([n n], ra2, rb3)
    sub2ind([n n], ra3, rb1)
    sub2ind([n n], ra3, rb2)
    sub2ind([n n], ra3, rb3)

    % Symmetric block: rows b, cols a
    sub2ind([n n], rb1, ra1)
    sub2ind([n n], rb1, ra2)
    sub2ind([n n], rb1, ra3)
    sub2ind([n n], rb2, ra1)
    sub2ind([n n], rb2, ra2)
    sub2ind([n n], rb2, ra3)
    sub2ind([n n], rb3, ra1)
    sub2ind([n n], rb3, ra2)
    sub2ind([n n], rb3, ra3)
];

val = [
    Txx
    Txy
    Txz
    Txy
    Tyy
    Tyz
    Txz
    Tyz
    Tzz

    Txx
    Txy
    Txz
    Txy
    Tyy
    Tyz
    Txz
    Tyz
    Tzz
];

Tpol(idx) = val;

nKept = numel(a);
fillTime = toc(tFill);

info = struct();
info.nPairBlocksKept = nKept;
info.nPairBlocksSkippedCutoff = nPol * (nPol - 1) / 2 - nKept;
info.cache_time = cacheTime;
info.fill_time = fillTime;
info.cache_info = local_summarize_cache(cache);

end

function [Tpol, info] = local_fill_all_pairs(Tpol, sys, activeSites, tensorOpts)

nPol = numel(activeSites);

nKept = 0;

for a = 1:(nPol - 1)
    i = activeSites(a);

    Ia = local_block_indices(a);

    for b = (a + 1):nPol
        j = activeSites(b);

        Ib = local_block_indices(b);

        Tij = thole.dipole_tensor_block( ...
            sys.site_pos(i, :), ...
            sys.site_pos(j, :), ...
            sys.site_alpha(i), ...
            sys.site_alpha(j), ...
            sys.thole_a, ...
            tensorOpts);

        Tpol(Ia, Ib) = Tij;
        Tpol(Ib, Ia) = Tij.';

        nKept = nKept + 1;
    end
end

info = struct();
info.nPairBlocksKept = nKept;
info.nPairBlocksSkippedCutoff = 0;
info.fill_time = NaN;

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

if ~isfield(sys, 'n_sites') || isempty(sys.n_sites)
    sys.n_sites = nSites; %#ok<NASGU>
end

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

function validate_options(useThole, softening, rcut, useMex, profile, verbose)

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

if ~(islogical(useMex) && isscalar(useMex))
    error('thole:assemble_nonperiodic_interaction_matrix:BadUseMex', ...
        'scfParams.use_mex must be a logical scalar.');
end

if ~(islogical(profile) && isscalar(profile))
    error('thole:assemble_nonperiodic_interaction_matrix:BadProfile', ...
        'scfParams.profile must be a logical scalar.');
end

if ~(islogical(verbose) && isscalar(verbose))
    error('thole:assemble_nonperiodic_interaction_matrix:BadVerbose', ...
        'scfParams.verbose must be a logical scalar.');
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

function summary = local_summarize_cache(cache)

summary = struct();

fields = {
    'n_pairs'
    'nPairs'
    'rcut'
    'used_mex'
    'use_mex'
    'backend'
    'build_time'
    'timing'
};

for k = 1:numel(fields)
    name = fields{k};

    if isfield(cache, name)
        summary.(name) = cache.(name);
    end
end

summary.n_pairs_returned = numel(cache.pair_i);

end