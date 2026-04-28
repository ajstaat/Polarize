function op = build_nonperiodic_paircache_operator(sys, problem, opt)
%BUILD_NONPERIODIC_PAIRCACHE_OPERATOR Build matrix-free pair-cache operator.
%
% Private backend builder used by thole.make_polarization_operator.
%
% Returns:
%   op.kind    = 'matrix_free'
%   op.backend = 'nonperiodic_paircache_apply'
%
% This backend stores pair-cache tensor coefficients and exposes:
%   op.apply(muVec)
%
% It does not build or store dense op.Tpol.

if ~isfinite(opt.Rcut)
    error('thole:build_nonperiodic_paircache_operator:RequiresFiniteCutoff', ...
        ['Nonperiodic pair-cache matrix-free operator requires finite Rcut. ', ...
         'Use Backend="dense" for all-pairs dense calculations.']);
end

io.assert_atomic_units(sys);

sites = problem.activeSites(:);
nPol = problem.nPolSites;
nSites = size(sys.site_pos, 1);

siteMask = false(nSites, 1);
siteMask(sites) = true;

cacheOpts = struct();
cacheOpts.rcut = opt.Rcut;
cacheOpts.site_mask = siteMask;
cacheOpts.use_mex = opt.UseMex;
cacheOpts.profile = opt.Profile || opt.Verbose;

tCache = tic;
cache = geom.build_nonperiodic_pair_cache(sys, cacheOpts);
cacheTime = toc(tCache);

fullToActive = zeros(nSites, 1);
fullToActive(sites) = 1:nPol;

pair_i_all = cache.pair_i(:);
pair_j_all = cache.pair_j(:);

a_all = fullToActive(pair_i_all);
b_all = fullToActive(pair_j_all);

if any(a_all == 0) || any(b_all == 0)
    error('thole:build_nonperiodic_paircache_operator:BadPairCacheMask', ...
        'Pair cache returned a site outside the active-site mask.');
end

keep = (a_all ~= b_all);

pair_i = pair_i_all(keep);
pair_j = pair_j_all(keep);
a = a_all(keep);
b = b_all(keep);

dr = cache.dr(keep, :);
r_bare = cache.r_bare(keep);
r2_bare = cache.r2_bare(keep);
inv_r3_bare = cache.inv_r3_bare(keep);
inv_r5_bare = cache.inv_r5_bare(keep);

dx = dr(:, 1);
dy = dr(:, 2);
dz = dr(:, 3);

if opt.Softening == 0
    invR3 = inv_r3_bare(:);
    invR5 = inv_r5_bare(:);
else
    r2 = r2_bare(:) + opt.Softening^2;
    r = sqrt(r2);
    invR3 = 1 ./ (r2 .* r);
    invR5 = invR3 ./ r2;
end

if opt.UseThole
    if isfield(cache, 'thole_f3') && isfield(cache, 'thole_f5')
        f3all = cache.thole_f3(:);
        f5all = cache.thole_f5(:);

        f3 = f3all(keep);
        f5 = f5all(keep);
    else
        tf = thole.thole_f3f5_factors( ...
            r_bare(:), ...
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

Txx = c5 .* dx .* dx - c3;
Txy = c5 .* dx .* dy;
Txz = c5 .* dx .* dz;

Tyy = c5 .* dy .* dy - c3;
Tyz = c5 .* dy .* dz;

Tzz = c5 .* dz .* dz - c3;

op = struct();
op.mode = 'nonperiodic';
op.kind = 'matrix_free';
op.backend = 'nonperiodic_paircache_apply';

op.nPolSites = nPol;
op.size = [3*nPol 3*nPol];

op.apply = @(muVec) local_apply_nonperiodic_paircache( ...
    muVec, nPol, a, b, Txx, Txy, Txz, Tyy, Tyz, Tzz);

op.cache = struct();
op.cache.source = 'geom.build_nonperiodic_pair_cache';
op.cache.pair_i = pair_i;
op.cache.pair_j = pair_j;
op.cache.active_i = a;
op.cache.active_j = b;
op.cache.n_pairs = numel(a);

op.info = struct();
op.info.assembly_backend = 'nonperiodic_paircache_apply';
op.info.nPolSites = nPol;
op.info.nPairBlocks = nPol * (nPol - 1) / 2;
op.info.nPairBlocksKept = numel(a);
op.info.nPairBlocksSkippedCutoff = op.info.nPairBlocks - numel(a);
op.info.rcut = opt.Rcut;
op.info.use_cutoff = true;
op.info.use_mex = opt.UseMex;
op.info.cache_time = cacheTime;
op.info.cache_info = local_summarize_cache(cache);

op.capabilities = struct();
op.capabilities.apply = true;
op.capabilities.dense_matrix = false;
op.capabilities.row_update = false;

op.params = struct();
op.params.use_thole = opt.UseThole;
op.params.softening = opt.Softening;
op.params.rcut = opt.Rcut;
op.params.use_mex = opt.UseMex;

end

% =========================================================================
% Local apply
% =========================================================================

function Tmu = local_apply_nonperiodic_paircache(muVec, nPol, a, b, Txx, Txy, Txz, Tyy, Tyz, Tzz)

muVec = muVec(:);

if numel(muVec) ~= 3*nPol
    error('thole:build_nonperiodic_paircache_operator:BadApplyVectorSize', ...
        'muVec must have length 3*nPolSites.');
end

mux = muVec(1:3:end);
muy = muVec(2:3:end);
muz = muVec(3:3:end);

% Contributions to active site a from active site b.
Eax = Txx .* mux(b) + Txy .* muy(b) + Txz .* muz(b);
Eay = Txy .* mux(b) + Tyy .* muy(b) + Tyz .* muz(b);
Eaz = Txz .* mux(b) + Tyz .* muy(b) + Tzz .* muz(b);

% Contributions to active site b from active site a.
Ebx = Txx .* mux(a) + Txy .* muy(a) + Txz .* muz(a);
Eby = Txy .* mux(a) + Tyy .* muy(a) + Tyz .* muz(a);
Ebz = Txz .* mux(a) + Tyz .* muy(a) + Tzz .* muz(a);

Ex = accumarray([a; b], [Eax; Ebx], [nPol 1], @sum, 0);
Ey = accumarray([a; b], [Eay; Eby], [nPol 1], @sum, 0);
Ez = accumarray([a; b], [Eaz; Ebz], [nPol 1], @sum, 0);

Tmu = zeros(3*nPol, 1);
Tmu(1:3:end) = Ex;
Tmu(2:3:end) = Ey;
Tmu(3:3:end) = Ez;

end

function summary = local_summarize_cache(cache)

summary = struct();

if isfield(cache, 'pair_i')
    summary.n_pairs_returned = numel(cache.pair_i);
end

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

end