function Edip = induced_field_from_dipoles_thole_periodic(sys, mu, ewaldParams, dipoleParams)
%INDUCED_FIELD_FROM_DIPOLES_THOLE_PERIODIC Periodic induced field from dipoles.
%
% Edip = thole.induced_field_from_dipoles_thole_periodic(sys, mu, ewaldParams)
% Edip = thole.induced_field_from_dipoles_thole_periodic(sys, mu, ewaldParams, dipoleParams)
%
% Inputs
%   sys           canonical polarization-system struct in atomic units
%   mu            N x 3 induced dipoles
%   ewaldParams   struct with fields:
%                   .alpha
%                   .rcut
%                   .kcut
%                   .boundary   optional, default 'tinfoil'
%
%   dipoleParams  optional struct with fields:
%                   .target_mask            N x 1 logical, optional
%                   .source_mask            N x 1 logical, optional
%                   .use_thole              logical, default true
%                   .realspace_cache        optional cache from
%                                           geom.build_periodic_realspace_cache(...)
%                   .kspace_cache           optional cache from
%                                           geom.build_periodic_kspace_cache(...)
%                   .problem                optional problem struct from
%                                           thole.prepare_scf_problem(...)
%
% Output
%   Edip          N x 3 periodic induced field
%
% Notes
%   The periodic operator is applied as:
%     E = E_real + E_recip + E_self + E_surf
%
%   Real-space action is cache-based with scalar coefficients:
%     E_{i<-j} = coeff_iso * mu_j + coeff_dyad * dr * (dr · mu_j)
%
%   Reciprocal-space action uses cached phases and k-space invariants:
%     E_i^recip = sum_k 2*pref(k) * [ cos_i * A_k + sin_i * B_k ] * (k k^T)
%   where
%     A_k = sum_j cos_j * (k·mu_j)
%     B_k = sum_j sin_j * (k·mu_j)
%
%   The analytic self term contributes
%     E_self(i) = -(4*alpha^3/(3*sqrt(pi))) * mu_i
%
%   For vacuum boundary conditions, the surface term contributes
%     E_surf(i) = (4*pi/(3V)) * sum_j mu_j
%   and is zero for tinfoil boundaries.

    io.assert_atomic_units(sys);

    if nargin < 4 || isempty(dipoleParams)
        dipoleParams = struct();
    end

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('thole:induced_field_from_dipoles_thole_periodic:MissingSitePos', ...
            'sys.site_pos is missing or empty.');
    end
    if ~isfield(sys, 'site_is_polarizable') || isempty(sys.site_is_polarizable)
        error('thole:induced_field_from_dipoles_thole_periodic:MissingSiteMask', ...
            'sys.site_is_polarizable is missing or empty.');
    end
    if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
        error('thole:induced_field_from_dipoles_thole_periodic:MissingSiteAlpha', ...
            'sys.site_alpha is missing or empty.');
    end

    nSites = size(sys.site_pos, 1);
    if ~isequal(size(mu), [nSites, 3])
        error('thole:induced_field_from_dipoles_thole_periodic:BadMuSize', ...
            'mu must be N x 3.');
    end

    if ~isfield(ewaldParams, 'alpha') || isempty(ewaldParams.alpha)
        error('thole:induced_field_from_dipoles_thole_periodic:MissingAlpha', ...
            'ewaldParams.alpha is required.');
    end
    if ~isfield(ewaldParams, 'rcut') || isempty(ewaldParams.rcut)
        error('thole:induced_field_from_dipoles_thole_periodic:MissingRcut', ...
            'ewaldParams.rcut is required.');
    end
    if ~isfield(ewaldParams, 'kcut') || isempty(ewaldParams.kcut)
        error('thole:induced_field_from_dipoles_thole_periodic:MissingKcut', ...
            'ewaldParams.kcut is required.');
    end

    target_mask = true(nSites, 1);
    if isfield(dipoleParams, 'target_mask') && ~isempty(dipoleParams.target_mask)
        target_mask = logical(dipoleParams.target_mask(:));
    end
    if numel(target_mask) ~= nSites
        error('thole:induced_field_from_dipoles_thole_periodic:BadTargetMask', ...
            'target_mask must have length N.');
    end

    source_mask = logical(sys.site_is_polarizable(:));
    if isfield(dipoleParams, 'source_mask') && ~isempty(dipoleParams.source_mask)
        source_mask = logical(dipoleParams.source_mask(:));
    end
    if numel(source_mask) ~= nSites
        error('thole:induced_field_from_dipoles_thole_periodic:BadSourceMask', ...
            'source_mask must have length N.');
    end

    use_thole = true;
    if isfield(dipoleParams, 'use_thole') && ~isempty(dipoleParams.use_thole)
        use_thole = logical(dipoleParams.use_thole);
    end

    problem = [];
    if isfield(dipoleParams, 'problem') && ~isempty(dipoleParams.problem)
        problem = dipoleParams.problem;
    end

    if isempty(problem)
        % Minimal problem-like struct for cache builders.
        polMask = logical(sys.site_is_polarizable(:));
        activeSites = find(polMask);
        problem = struct();
        problem.nSites = nSites;
        problem.polMask = polMask;
        problem.activeSites = activeSites;
        problem.nPolSites = numel(activeSites);
        problem.activeVecIdx = zeros(3 * numel(activeSites), 1); %#ok<STRNU>
    end

    realCache = [];
    if isfield(dipoleParams, 'realspace_cache') && ~isempty(dipoleParams.realspace_cache)
        realCache = dipoleParams.realspace_cache;
    else
        cacheParams = struct();
        cacheParams.use_thole = use_thole;
        realCache = geom.build_periodic_realspace_cache(sys, problem, ewaldParams, cacheParams);
    end

    kCache = [];
    if isfield(dipoleParams, 'kspace_cache') && ~isempty(dipoleParams.kspace_cache)
        kCache = dipoleParams.kspace_cache;
    else
        kCache = geom.build_periodic_kspace_cache(sys, problem, ewaldParams);
    end

    Edip = zeros(nSites, 3);

    % ---------------------------------------------------------------------
    % Real-space contribution

    pair_i = realCache.pair_i(:);
    pair_j = realCache.pair_j(:);
    dr = realCache.dr;
    coeff_iso = realCache.coeff_iso(:);
    coeff_dyad = realCache.coeff_dyad(:);

    nReal = realCache.nInteractions;

    for p = 1:nReal
        i = pair_i(p);
        j = pair_j(p);

        if i == j
            % Self-image interaction: one directed contribution on-site.
            if ~(target_mask(i) && source_mask(i))
                continue;
            end

            muj = mu(i, :);
            if all(muj == 0)
                continue;
            end

            rij = dr(p, :);   % image displacement
            muDotR = dot(muj, rij);

            Edip(i, :) = Edip(i, :) ...
                + coeff_iso(p) * muj ...
                + coeff_dyad(p) * rij * muDotR;
        else
            rij = dr(p, :);   % r_j + R - r_i

            % Field at i from source j
            if target_mask(i) && source_mask(j)
                muj = mu(j, :);
                if ~all(muj == 0)
                    muDotR = dot(muj, rij);
                    Edip(i, :) = Edip(i, :) ...
                        + coeff_iso(p) * muj ...
                        + coeff_dyad(p) * rij * muDotR;
                end
            end

            % Field at j from source i, using reversed directed displacement
            if target_mask(j) && source_mask(i)
                mui = mu(i, :);
                if ~all(mui == 0)
                    rji = -rij;
                    muDotR = dot(mui, rji);
                    Edip(j, :) = Edip(j, :) ...
                        + coeff_iso(p) * mui ...
                        + coeff_dyad(p) * rji * muDotR;
                end
            end
        end
    end

    % ---------------------------------------------------------------------
    % Reciprocal-space contribution
    %
    % For each k:
    %   v_j = k · mu_j
    %   A_k = sum_j cos_j * v_j
    %   B_k = sum_j sin_j * v_j
    %
    % Then for each target i:
    %   phase_factor_i = cos_i * A_k + sin_i * B_k
    %   E_i += 2 * pref(k) * phase_factor_i * (k k^T)
    %
    % Here we evaluate targets in active-space form, then scatter back.

    activeSites = kCache.activeSites(:);
    nPol = kCache.nPolSites;

    if kCache.num_kvec > 0 && nPol > 0
        mu_pol = mu(activeSites, :);

        % Apply source mask in active space.
        source_active = source_mask(activeSites);
        mu_src = mu_pol;
        mu_src(~source_active, :) = 0;

        cos_phase = kCache.cos_phase;   % nPol x Nk
        sin_phase = kCache.sin_phase;   % nPol x Nk
        kvecs = kCache.kvecs;           % Nk x 3
        pref = kCache.pref(:);          % Nk x 1
        kk6 = kCache.kk6;               %#ok<NASGU>

        % v = k · mu_j for all active sites and all k
        % mu_src: nPol x 3, kvecs': 3 x Nk => nPol x Nk
        v = mu_src * kvecs.';

        % Structure-factor-like sums over sources
        A = sum(cos_phase .* v, 1);   % 1 x Nk
        B = sum(sin_phase .* v, 1);   % 1 x Nk

        % Phase factor for each target site and each k
        phase_factor = cos_phase .* A + sin_phase .* B;   % nPol x Nk

        % Weight each k contribution
        W = phase_factor .* (2 * pref.');                 % nPol x Nk

        % Target-wise reciprocal field components:
        %   Ex = sum_k W .* kx .* (k·e_x? no) => W * [kxx, kxy, kxz]
        % More directly:
        %   E_i = sum_k W(i,k) * [kxx kxy kxz; kxy kyy kyz; kxz kyz kzz] * [1;1;1]? no
        % Since reciprocal operator acting on mu is already condensed into W through
        % the source contraction v = k·mu, the field is:
        %   E_i = sum_k W(i,k) * kvec(k,:)
        % Wait: T_ij = 2*pref*cos * (k k^T), so acting on mu_j:
        %   (k k^T) mu_j = k * (k·mu_j)
        % Summing over j created A/B from (k·mu_j), so the target field is
        %   E_i = sum_k 2*pref*phase_factor_i * k
        % component-wise multiplied by k.
        %
        % Thus:
        Ex_pol = W * kvecs(:, 1);
        Ey_pol = W * kvecs(:, 2);
        Ez_pol = W * kvecs(:, 3);

        Erecip_pol = [Ex_pol, Ey_pol, Ez_pol];

        % Apply target mask in active space before scatter.
        target_active = target_mask(activeSites);
        Erecip_pol(~target_active, :) = 0;

        Edip(activeSites, :) = Edip(activeSites, :) + Erecip_pol;
    end

    % ---------------------------------------------------------------------
    % Analytic self term
    %
    % Tself block = -(4*alpha^3/(3*sqrt(pi))) I
    % so E_self(i) = Tself * mu_i

    alpha = ewaldParams.alpha;
    self_coeff = -(4 * alpha^3 / (3 * sqrt(pi)));

    self_sites = target_mask & source_mask;
    Edip(self_sites, :) = Edip(self_sites, :) + self_coeff * mu(self_sites, :);

    % ---------------------------------------------------------------------
    % Surface term
    %
    % For vacuum:
    %   Tsurf block = (4*pi/(3V)) I between every pair
    % so
    %   E_surf(i) = (4*pi/(3V)) * sum_j mu_j
    % over allowed source sites.
    %
    % For tinfoil: zero.

    boundary = 'tinfoil';
    if isfield(ewaldParams, 'boundary') && ~isempty(ewaldParams.boundary)
        boundary = lower(char(string(ewaldParams.boundary)));
    end

    switch boundary
        case 'tinfoil'
            % no-op
        case 'vacuum'
            H = local_get_direct_lattice(sys);
            V = abs(det(H));
            surf_coeff = 4 * pi / (3 * V);

            Msrc = sum(mu(source_mask, :), 1);
            Edip(target_mask, :) = Edip(target_mask, :) + surf_coeff * Msrc;
        otherwise
            error('thole:induced_field_from_dipoles_thole_periodic:UnknownBoundary', ...
                'boundary must be ''tinfoil'' or ''vacuum''.');
    end
end

function H = local_get_direct_lattice(sys)
%LOCAL_GET_DIRECT_LATTICE Extract direct lattice matrix, columns are lattice vectors.

    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('thole:induced_field_from_dipoles_thole_periodic:MissingLattice', ...
            'Need sys.super_lattice or sys.lattice.');
    end

    validateattributes(H, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
        mfilename, 'lattice');
end