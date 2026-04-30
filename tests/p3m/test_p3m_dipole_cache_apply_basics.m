function test_p3m_dipole_cache_apply_basics()
%TEST_P3M_DIPOLE_CACHE_APPLY_BASICS Basic P3M dipole-cache tests.
%
% Covers:
%   p3m.build_dipole_cache
%   p3m.apply_dipole_cache
%
% Checks:
%   - cache builds with project row-lattice convention
%   - zero dipoles produce exactly zero reciprocal field
%   - full-system N x 3, source-local Nsrc x 3, and stacked 3*Nsrc inputs agree
%   - assigned mesh dipole moment is conserved in cache apply diagnostics
%   - integer-cell shifted input geometry gives the same cached apply result
%
% This tests the production induced-dipole P3M route, not a one-shot
% duplicate induced-dipole implementation.

rng(41);

sys = local_make_periodic_polsys();

opts = struct();
opts.ewald = struct();
opts.ewald.alpha = 0.30;
opts.mesh_size = [12 10 8];
opts.assignment_order = 4;
opts.target_mask = logical(sys.site_is_polarizable(:));
opts.source_mask = logical(sys.site_is_polarizable(:));
opts.derivative_mode = 'spectral';
opts.influence_mode = 'ewald';
opts.deconvolve_assignment = true;
opts.deconvolution_floor = 1e-8;
opts.verbose = false;

cache = p3m.build_dipole_cache(sys, opts);

assert(strcmp(cache.mode, 'p3m_dipole_reciprocal_cache'), ...
    'Unexpected P3M dipole cache mode.');
assert(strcmp(cache.lattice_convention, 'project_row_H_column_G_HG_2piI'), ...
    'P3M dipole cache should use canonical lattice convention.');
assert(norm(cache.H * cache.G - 2*pi*eye(3), 'fro') < 1e-12, ...
    'P3M dipole cache should satisfy H*G = 2*pi*I.');
assert(cache.nSources == nnz(opts.source_mask), ...
    'Cache nSources should match source mask.');
assert(cache.nTargets == nnz(opts.target_mask), ...
    'Cache nTargets should match target mask.');
assert(cache.nK > 0, ...
    'P3M dipole cache should include nonzero reciprocal modes.');

nSites = sys.n_sites;
nSources = cache.nSources;

%% ------------------------------------------------------------------------
% Zero dipoles
% -------------------------------------------------------------------------

muZero = zeros(nSites, 3);

[Ezero, partsZero] = p3m.apply_dipole_cache(cache, muZero);

assert(isequal(size(Ezero), [nSites, 3]), ...
    'P3M cached apply should return Nsites x 3 field.');
assert(norm(Ezero(:)) == 0, ...
    'Zero dipoles should produce exactly zero cached P3M reciprocal field.');
assert(norm(partsZero.P_total) == 0, ...
    'Zero dipoles should assign zero total mesh dipole.');

%% ------------------------------------------------------------------------
% Equivalent input forms
% -------------------------------------------------------------------------

muFull = zeros(nSites, 3);
muSource = [
     0.10  -0.20   0.30
    -0.40   0.50  -0.60
     0.70   0.10  -0.20
    -0.30  -0.10   0.40
];

muFull(cache.source_sites, :) = muSource;
muVec = util.stack_xyz(muSource);

[Efull, partsFull] = p3m.apply_dipole_cache(cache, muFull);
[Esource, partsSource] = p3m.apply_dipole_cache(cache, muSource);
[Evec, partsVec] = p3m.apply_dipole_cache(cache, muVec);

assert(all(isfinite(Efull(:))), ...
    'Cached P3M reciprocal field should be finite.');
assert(norm(Efull(opts.target_mask, :), 'fro') > 0, ...
    'Nonzero source dipoles should produce nonzero target field.');

assert(norm(Efull - Esource, 'fro') < 1e-13 * max(1, norm(Efull, 'fro')), ...
    'Full-system and source-local matrix dipole inputs should agree.');
assert(norm(Efull - Evec, 'fro') < 1e-13 * max(1, norm(Efull, 'fro')), ...
    'Full-system and stacked-vector dipole inputs should agree.');

assert(norm(partsFull.P_total - sum(muSource, 1)) < 1e-13, ...
    'Cached scatter should conserve total source dipole.');
assert(norm(partsSource.P_total - sum(muSource, 1)) < 1e-13, ...
    'Source-local cached scatter should conserve total source dipole.');
assert(norm(partsVec.P_total - sum(muSource, 1)) < 1e-13, ...
    'Stacked-vector cached scatter should conserve total source dipole.');

assert(norm(partsFull.P_total - partsSource.P_total) < 1e-13, ...
    'Equivalent inputs should report identical total mesh dipoles.');
assert(norm(partsFull.P_total - partsVec.P_total) < 1e-13, ...
    'Equivalent inputs should report identical total mesh dipoles.');

%% ------------------------------------------------------------------------
% Target masking
% -------------------------------------------------------------------------

nonTarget = ~opts.target_mask;

assert(norm(Efull(nonTarget, :), 'fro') == 0, ...
    'Cached P3M reciprocal field should be zero on non-target rows.');

%% ------------------------------------------------------------------------
% Integer-cell shift invariance
% -------------------------------------------------------------------------

sysShift = sys;
lat = geom.get_lattice(sys);

% Shift every site by an integer lattice vector. Fractional coordinates
% change by integers, so P3M assignment/interpolation should be identical.
integerShift = [1 -2 3];
cartShift = integerShift * lat.H;

sysShift.site_pos = sys.site_pos + cartShift;

cacheShift = p3m.build_dipole_cache(sysShift, opts);
[Eshift, partsShift] = p3m.apply_dipole_cache(cacheShift, muFull);

assert(norm(Efull - Eshift, 'fro') < 1e-12 * max(1, norm(Efull, 'fro')), ...
    'Cached P3M apply should be invariant to integer-cell translations.');

assert(norm(partsShift.P_total - partsFull.P_total) < 1e-13, ...
    'Integer-cell shifted cache should conserve the same total mesh dipole.');

io.assert_atomic_units(sys);
end

function sys = local_make_periodic_polsys()
sys = struct();

% Four polarizable source/target sites plus one spectator.
sys.site_pos = [
     3.0   2.0   2.0
    10.0   6.0   5.0
    17.0  11.0   9.0
     8.0  19.0  14.0
     6.0  15.0  12.0
];

sys.site_charge = zeros(5, 1);

sys.site_alpha = [
    0.08
    0.07
    0.06
    0.05
    0.00
];

sys.site_is_polarizable = [
    true
    true
    true
    true
    false
];

sys.n_sites = 5;
sys.thole_a = 0.39;

sys.site_type = {'X'; 'X'; 'X'; 'X'; 'X'};
sys.site_class = {'pol'; 'pol'; 'pol'; 'pol'; 'spectator'};
sys.site_label = {'p1'; 'p2'; 'p3'; 'p4'; 'n'};
sys.site_mol_id = [1; 2; 3; 4; 5];
sys.site_is_active = sys.site_is_polarizable;

sys.units = struct();
sys.units.length = 'bohr';
sys.units.alpha = 'atomic_unit';
sys.units.charge = 'elementary_charge';

sys.is_periodic = true;
sys.periodic_mode = 'periodic';

% Row-lattice convention: cart = frac * H.
sys.lattice = [
    30.0   0.0   0.0
     2.0  29.0   0.0
     1.0   3.0  28.0
];

sys.super_lattice = sys.lattice;
end