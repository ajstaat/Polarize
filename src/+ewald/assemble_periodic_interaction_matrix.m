function [Tpol, parts, opinfo] = assemble_periodic_interaction_matrix(sys, problem, ewaldParams, scfParams)
%ASSEMBLE_PERIODIC_INTERACTION_MATRIX Assemble triclinic periodic dipole operator
% in polarizable-only active space using split real-space and k-space caches.
%
% Inputs
%   sys         canonical polarization-system struct with fields:
%                 .site_pos
%                 .site_alpha
%                 .site_is_polarizable
%                 .super_lattice or .lattice
%
%   problem     struct from thole.prepare_scf_problem(...)
%
%   ewaldParams struct with fields:
%                 .alpha
%                 .rcut
%                 .kcut
%                 .boundary              optional, default 'tinfoil'
%
%   scfParams   optional struct with fields:
%                 .use_thole             logical, default true
%                 .verbose               logical, default false
%                 .printEvery            positive integer, default 25
%
% Outputs
%   Tpol        3Np x 3Np periodic interaction matrix in polarizable-only space
%   parts       struct with component matrices:
%                 .real
%                 .recip
%                 .self
%                 .surf
%   opinfo      metadata describing the active-space operator

    io.assert_atomic_units(sys);

    if nargin < 4 || isempty(scfParams)
        scfParams = struct();
    end

    if nargin < 3 || isempty(ewaldParams)
        error('ewald:assemble_periodic_interaction_matrix:MissingEwaldParams', ...
            'ewaldParams with alpha, rcut, and kcut is required.');
    end
    if nargin < 2 || isempty(problem)
        error('ewald:assemble_periodic_interaction_matrix:MissingProblem', ...
            'problem struct from thole.prepare_scf_problem is required.');
    end

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('ewald:assemble_periodic_interaction_matrix:MissingSitePos', ...
            'sys.site_pos is missing or empty.');
    end
    if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
        error('ewald:assemble_periodic_interaction_matrix:MissingSiteAlpha', ...
            'sys.site_alpha is missing or empty.');
    end
    if ~isfield(sys, 'site_is_polarizable') || isempty(sys.site_is_polarizable)
        error('ewald:assemble_periodic_interaction_matrix:MissingSiteMask', ...
            'sys.site_is_polarizable is missing or empty.');
    end
    if ~isfield(problem, 'activeSites') || ~isfield(problem, 'nPolSites')
        error('ewald:assemble_periodic_interaction_matrix:BadProblem', ...
            'problem must include activeSites and nPolSites.');
    end

    if ~isfield(ewaldParams, 'alpha') || isempty(ewaldParams.alpha)
        error('ewald:assemble_periodic_interaction_matrix:MissingAlpha', ...
            'ewaldParams.alpha is required.');
    end
    if ~isfield(ewaldParams, 'rcut') || isempty(ewaldParams.rcut)
        error('ewald:assemble_periodic_interaction_matrix:MissingRcut', ...
            'ewaldParams.rcut is required.');
    end
    if ~isfield(ewaldParams, 'kcut') || isempty(ewaldParams.kcut)
        error('ewald:assemble_periodic_interaction_matrix:MissingKcut', ...
            'ewaldParams.kcut is required.');
    end

    use_thole = true;
    if isfield(scfParams, 'use_thole') && ~isempty(scfParams.use_thole)
        use_thole = logical(scfParams.use_thole);
    end

    verbose = false;
    if isfield(scfParams, 'verbose') && ~isempty(scfParams.verbose)
        verbose = logical(scfParams.verbose);
    end

    printEvery = 25;
    if isfield(scfParams, 'printEvery') && ~isempty(scfParams.printEvery)
        printEvery = scfParams.printEvery;
    end

    sites = problem.activeSites(:);
    nPol = problem.nPolSites;
    nSites = problem.nSites;

    H = local_get_direct_lattice(sys);
    boundary = local_get_boundary(ewaldParams);

    if verbose
        fprintf('Periodic interaction matrix: %d total sites, %d polarizable sites\n', ...
            nSites, nPol);
        fprintf('Compressed matrix size: %d x %d\n', 3*nPol, 3*nPol);
        fprintf('Ewald parameters: alpha=%.6g, rcut=%.6g, kcut=%.6g, boundary=%s\n', ...
            ewaldParams.alpha, ewaldParams.rcut, ewaldParams.kcut, boundary);
    end

    % ---------------------------------------------------------------------
    % Build split caches

    realGeom = geom.build_periodic_realspace_geom_cache(sys, problem, ewaldParams);
    realRow  = geom.build_periodic_realspace_row_geom_cache(sys, problem, realGeom);
    realCoeff = geom.build_periodic_realspace_dipole_coeff_cache(sys, realRow, ewaldParams, scfParams);

    targetCache = geom.build_periodic_kspace_target_cache(sys, problem, ewaldParams);
    kCoeff = geom.build_periodic_kspace_dipole_coeff_cache(targetCache);

    % ---------------------------------------------------------------------
    % Allocate component matrices

    Treal  = zeros(3*nPol, 3*nPol);
    Trecip = zeros(3*nPol, 3*nPol);
    Tself  = zeros(3*nPol, 3*nPol);
    Tsurf  = zeros(3*nPol, 3*nPol);

    % ---------------------------------------------------------------------
    % Real-space contribution from row cache

    tReal = tic;

    row_ptr = realRow.row_ptr;
    col_idx = realRow.col_idx;
    dr = realRow.dr;
    coeff_iso = realCoeff.coeff_iso(:);
    coeff_dyad = realCoeff.coeff_dyad(:);

    I3 = eye(3);

    for a = 1:nPol
        Ia = util.block3(a);
        idx0 = row_ptr(a);
        idx1 = row_ptr(a + 1) - 1;

        if idx1 < idx0
            continue;
        end

        for p = idx0:idx1
            b = col_idx(p);
            Ib = util.block3(b);

            x = dr(p, :).';
            xxT = x * x.';
            Tij = coeff_iso(p) * I3 + coeff_dyad(p) * xxT;

            Treal(Ia, Ib) = Treal(Ia, Ib) + Tij;
        end

        if verbose && (a == 1 || a == nPol || mod(a, printEvery) == 0)
            frac      = a / max(nPol, 1);
            elapsed   = toc(tReal);
            totalEst  = elapsed / max(frac, eps);
            remainEst = max(totalEst - elapsed, 0);
            fprintf(' assembling Treal: row %d / %d (%.1f%%) elapsed %.1f s ETA %.1f s\n', ...
                a, nPol, 100*frac, elapsed, remainEst);
        end
    end

    % ---------------------------------------------------------------------
    % Reciprocal-space contribution from split target cache + dipole coeff cache

    tRecip = tic;

    if targetCache.num_kvec > 0
        nk = targetCache.num_kvec;
        pref = kCoeff.pref(:);
        kk6 = kCoeff.kk6;
        cos_phase = targetCache.target_cos;
        sin_phase = targetCache.target_sin;

        Txx = zeros(nPol, nPol);
        Tyy = zeros(nPol, nPol);
        Tzz = zeros(nPol, nPol);
        Txy = zeros(nPol, nPol);
        Txz = zeros(nPol, nPol);
        Tyz = zeros(nPol, nPol);

        blockSize = 256;

        for k0 = 1:blockSize:nk
            k1 = min(k0 + blockSize - 1, nk);
            idx = k0:k1;

            C = cos_phase(:, idx);
            S = sin_phase(:, idx);

            Txx = Txx + local_accumulate_component(C, S, 2 * pref(idx) .* kk6(idx,1));
            Tyy = Tyy + local_accumulate_component(C, S, 2 * pref(idx) .* kk6(idx,2));
            Tzz = Tzz + local_accumulate_component(C, S, 2 * pref(idx) .* kk6(idx,3));
            Txy = Txy + local_accumulate_component(C, S, 2 * pref(idx) .* kk6(idx,4));
            Txz = Txz + local_accumulate_component(C, S, 2 * pref(idx) .* kk6(idx,5));
            Tyz = Tyz + local_accumulate_component(C, S, 2 * pref(idx) .* kk6(idx,6));

            if verbose && (k0 == 1 || k1 == nk || mod(k1, printEvery * blockSize) == 0)
                frac      = k1 / max(nk, 1);
                elapsed   = toc(tRecip);
                totalEst  = elapsed / max(frac, eps);
                remainEst = max(totalEst - elapsed, 0);
                fprintf(' assembling Trecip: k %5d / %5d (%.1f%%) elapsed %.1f s ETA %.1f s\n', ...
                    k1, nk, 100*frac, elapsed, remainEst);
            end
        end

        for a = 1:nPol
            Ia = util.block3(a);
            for b = 1:nPol
                Ib = util.block3(b);
                Trecip(Ia, Ib) = [ ...
                    Txx(a,b), Txy(a,b), Txz(a,b); ...
                    Txy(a,b), Tyy(a,b), Tyz(a,b); ...
                    Txz(a,b), Tyz(a,b), Tzz(a,b)];
            end
        end
    end

    % ---------------------------------------------------------------------
    % Analytic self term

    selfBlock = ewald.self_tensor_block_dipole(ewaldParams.alpha);
    for a = 1:nPol
        Ia = util.block3(a);
        Tself(Ia, Ia) = selfBlock;
    end

    % ---------------------------------------------------------------------
    % Surface term

    surfBlock = ewald.surface_tensor_block_dipole(H, boundary);
    if any(surfBlock(:) ~= 0)
        for a = 1:nPol
            Ia = util.block3(a);
            for b = 1:nPol
                Ib = util.block3(b);
                Tsurf(Ia, Ib) = surfBlock;
            end
        end
    end

    % ---------------------------------------------------------------------
    % Final operator and metadata

    Tpol = Treal + Trecip + Tself + Tsurf;

    parts = struct();
    parts.real  = Treal;
    parts.recip = Trecip;
    parts.self  = Tself;
    parts.surf  = Tsurf;

    opinfo = struct();
    opinfo.mode = 'periodic';
    opinfo.space = 'polarizable_only';
    opinfo.nSites = nSites;
    opinfo.nPolSites = nPol;
    opinfo.activeSites = sites;
    opinfo.activeVecIdx = problem.activeVecIdx;

    opinfo.alpha = ewaldParams.alpha;
    opinfo.rcut = ewaldParams.rcut;
    opinfo.kcut = ewaldParams.kcut;
    opinfo.boundary = boundary;
    opinfo.use_thole_real_space = use_thole;

    opinfo.real_image_bounds = realGeom.real_image_bounds;
    opinfo.nRealInteractions = realGeom.nInteractions;
    opinfo.num_kvec = targetCache.num_kvec;
    opinfo.volume = targetCache.V;

    if verbose
        fprintf('Finished assembling compressed periodic interaction matrix.\n');
        fprintf('  nRealInteractions = %d\n', realGeom.nInteractions);
        fprintf('  num_kvec          = %d\n', targetCache.num_kvec);
    end
end

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('ewald:assemble_periodic_interaction_matrix:MissingLattice', ...
            'Need sys.super_lattice or sys.lattice.');
    end

    validateattributes(H, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
        mfilename, 'lattice');
end

function boundary = local_get_boundary(ewaldParams)
    boundary = 'tinfoil';
    if isfield(ewaldParams, 'boundary') && ~isempty(ewaldParams.boundary)
        if ~(ischar(ewaldParams.boundary) || isstring(ewaldParams.boundary))
            error('ewald:assemble_periodic_interaction_matrix:BadBoundary', ...
                'ewaldParams.boundary must be a character vector or string scalar.');
        end
        boundary = lower(char(string(ewaldParams.boundary)));
    end
end

function T = local_accumulate_component(C, S, d)
    d = d(:).';

    pos = d > 0;
    neg = d < 0;

    T = zeros(size(C,1), size(C,1));

    if any(pos)
        wp = sqrt(d(pos));
        Cp = C(:, pos) .* wp;
        Sp = S(:, pos) .* wp;
        T = T + Cp * Cp.' + Sp * Sp.';
    end

    if any(neg)
        wn = sqrt(-d(neg));
        Cn = C(:, neg) .* wn;
        Sn = S(:, neg) .* wn;
        T = T - (Cn * Cn.' + Sn * Sn.');
    end
end