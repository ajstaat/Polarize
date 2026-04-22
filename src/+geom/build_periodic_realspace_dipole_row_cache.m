function rowCache = build_periodic_realspace_dipole_row_cache(sys, problem, ewaldParams, scfParams)
%BUILD_PERIODIC_REALSPACE_DIPOLE_ROW_CACHE Build direct solver-facing periodic
% real-space dipole row cache.
%
% rowCache = geom.build_periodic_realspace_dipole_row_cache(sys, problem, ewaldParams)
% rowCache = geom.build_periodic_realspace_dipole_row_cache(sys, problem, ewaldParams, scfParams)
%
% Fast path:
%   - uses mex_build_periodic_realspace_dipole_row_cache when available
%
% Fallback:
%   - preserves the existing layered periodic workflow:
%       geom cache -> row geom cache -> dipole coeff cache
%
% Output fields
%   .mode
%   .nSites
%   .nPolSites
%   .activeSites
%   .targetSites
%   .sourceSites
%   .full_to_active
%   .full_to_target
%   .full_to_source
%   .row_ptr
%   .col_idx
%   .source_full_idx
%   .dr
%   .r_bare
%   .r2_bare
%   .coeff_iso
%   .coeff_dyad
%   .alpha
%   .rcut
%   .use_thole
%   .nInteractions
%   .real_image_bounds
%   .nImageShifts
%   .timing   (if profile requested)

    narginchk(3, 4);

    if nargin < 4 || isempty(scfParams)
        scfParams = struct();
    end

    validateattributes(problem, {'struct'}, {'scalar'}, mfilename, 'problem', 2);
    validateattributes(ewaldParams, {'struct'}, {'scalar'}, mfilename, 'ewaldParams', 3);

    io.assert_atomic_units(sys);

    if ~isfield(problem, 'activeSites') || isempty(problem.activeSites)
        error('geom:build_periodic_realspace_dipole_row_cache:MissingActiveSites', ...
            'problem.activeSites is required.');
    end
    activeSites = problem.activeSites(:);
    nPolSites = numel(activeSites);
    nSites = problem.nSites;

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('geom:build_periodic_realspace_dipole_row_cache:MissingSitePos', ...
            'sys.site_pos is required.');
    end
    pos = sys.site_pos;
    validateattributes(pos, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
        mfilename, 'sys.site_pos');
    if size(pos, 1) ~= nSites
        error('geom:build_periodic_realspace_dipole_row_cache:BadSitePos', ...
            'size(sys.site_pos,1) must match problem.nSites.');
    end

    if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
        error('geom:build_periodic_realspace_dipole_row_cache:MissingSiteAlpha', ...
            'sys.site_alpha is required.');
    end
    alphaSite = sys.site_alpha(:);
    if numel(alphaSite) ~= nSites
        error('geom:build_periodic_realspace_dipole_row_cache:BadSiteAlpha', ...
            'sys.site_alpha must have length nSites.');
    end

    if ~isfield(ewaldParams, 'alpha') || isempty(ewaldParams.alpha)
        error('geom:build_periodic_realspace_dipole_row_cache:MissingAlpha', ...
            'ewaldParams.alpha is required.');
    end
    alpha = ewaldParams.alpha;
    validateattributes(alpha, {'double'}, {'scalar','real','finite','positive'}, ...
        mfilename, 'ewaldParams.alpha');

    if ~isfield(ewaldParams, 'rcut') || isempty(ewaldParams.rcut)
        error('geom:build_periodic_realspace_dipole_row_cache:MissingRcut', ...
            'ewaldParams.rcut is required.');
    end
    rcut = ewaldParams.rcut;
    validateattributes(rcut, {'double'}, {'scalar','real','finite','positive'}, ...
        mfilename, 'ewaldParams.rcut');

    use_thole = true;
    if isfield(scfParams, 'use_thole') && ~isempty(scfParams.use_thole)
        use_thole = logical(scfParams.use_thole);
    end

    thole_a = 0.0;
    if use_thole
        if ~isfield(sys, 'thole_a') || isempty(sys.thole_a)
            error('geom:build_periodic_realspace_dipole_row_cache:MissingTholeA', ...
                'sys.thole_a is required when use_thole = true.');
        end
        thole_a = sys.thole_a;
    end

    doProfile = false;
    if isfield(scfParams, 'profile') && ~isempty(scfParams.profile)
        doProfile = logical(scfParams.profile);
    end

    H = local_get_direct_lattice(sys);
    a = H(:,1);
    b = H(:,2);
    c = H(:,3);
    nxmax = ceil(rcut / norm(a)) + 1;
    nymax = ceil(rcut / norm(b)) + 1;
    nzmax = ceil(rcut / norm(c)) + 1;
    nImageShifts = (2*nxmax + 1) * (2*nymax + 1) * (2*nzmax + 1);

    mexAvailable = (exist(['mex_build_periodic_realspace_dipole_row_cache.' mexext], 'file') == 3) || ...
                   (exist('mex_build_periodic_realspace_dipole_row_cache', 'file') == 3);

    tTotal = tic;
    timing = struct();

    if mexAvailable
        t0 = tic;
        [row_ptr, col_idx, source_full_idx, dr, r_bare, r2_bare, coeff_iso, coeff_dyad, ...
            fullToActive, nDirected] = mex_build_periodic_realspace_dipole_row_cache( ...
                double(pos), ...
                double(activeSites), ...
                double(H), ...
                double(rcut), ...
                double(alpha), ...
                double(alphaSite), ...
                double(thole_a), ...
                double(use_thole));
        timing.mex_build_s = toc(t0);

        rowCache = struct();
        rowCache.mode = 'periodic_realspace_dipole_row';
        rowCache.nSites = nSites;
        rowCache.nPolSites = nPolSites;
        rowCache.activeSites = activeSites;
        rowCache.targetSites = activeSites;
        rowCache.sourceSites = activeSites;
        rowCache.full_to_active = fullToActive;
        rowCache.full_to_target = fullToActive;
        rowCache.full_to_source = fullToActive;

        rowCache.row_ptr = row_ptr;
        rowCache.col_idx = col_idx;
        rowCache.source_full_idx = source_full_idx;
        rowCache.dr = dr;
        rowCache.r_bare = r_bare;
        rowCache.r2_bare = r2_bare;
        rowCache.coeff_iso = coeff_iso;
        rowCache.coeff_dyad = coeff_dyad;

        rowCache.alpha = alpha;
        rowCache.rcut = rcut;
        rowCache.use_thole = use_thole;
        rowCache.nInteractions = nDirected;
        rowCache.real_image_bounds = [nxmax, nymax, nzmax];
        rowCache.nImageShifts = nImageShifts;

        timing.total_s = toc(tTotal);
        if doProfile
            rowCache.timing = timing;
            fprintf('build_periodic_realspace_dipole_row_cache timing summary:\n');
            fprintf('  mex_build            : %.6f s\n', timing.mex_build_s);
            fprintf('  total                : %.6f s\n', timing.total_s);
            fprintf('  nInteractions        : %d\n', rowCache.nInteractions);
            fprintf('  nImageShifts         : %d\n', rowCache.nImageShifts);
            fprintf('  backend              : MEX\n');
        end
        return;
    end

    % ------------------------------------------------------------------
    % MATLAB fallback via existing layered periodic caches

    geomOpts = struct();
    geomOpts.target_mask = problem.polMask;
    geomOpts.source_mask = problem.polMask;

    t0 = tic;
    geomCache = geom.build_periodic_realspace_geom_cache(sys, problem, ewaldParams, geomOpts);
    timing.geom_s = toc(t0);

    t0 = tic;
    rowGeomCache = geom.build_periodic_realspace_row_geom_cache(sys, problem, geomCache);
    timing.row_geom_s = toc(t0);

    t0 = tic;
    coeffCache = geom.build_periodic_realspace_dipole_coeff_cache(sys, rowGeomCache, ewaldParams, scfParams);
    timing.coeff_s = toc(t0);

    rowCache = rowGeomCache;
    rowCache.mode = 'periodic_realspace_dipole_row';
    rowCache.coeff_iso = coeffCache.coeff_iso;
    rowCache.coeff_dyad = coeffCache.coeff_dyad;
    rowCache.alpha = coeffCache.alpha;
    rowCache.use_thole = coeffCache.use_thole;
    rowCache.real_image_bounds = geomCache.real_image_bounds;
    rowCache.nImageShifts = geomCache.nImageShifts;

    timing.total_s = toc(tTotal);
    if doProfile
        rowCache.timing = timing;
        fprintf('build_periodic_realspace_dipole_row_cache timing summary:\n');
        fprintf('  geom_cache           : %.6f s\n', timing.geom_s);
        fprintf('  row_geom_cache       : %.6f s\n', timing.row_geom_s);
        fprintf('  dipole_coeff_cache   : %.6f s\n', timing.coeff_s);
        fprintf('  total                : %.6f s\n', timing.total_s);
        fprintf('  nInteractions        : %d\n', rowCache.nInteractions);
        fprintf('  nImageShifts         : %d\n', rowCache.nImageShifts);
        fprintf('  backend              : MATLAB layered fallback\n');
    end
end

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('geom:build_periodic_realspace_dipole_row_cache:MissingLattice', ...
            'Missing direct lattice on system.');
    end
end