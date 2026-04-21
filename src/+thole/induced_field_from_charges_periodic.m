function E = induced_field_from_charges_periodic(sys, fieldParams)
    if nargin < 2 || isempty(fieldParams)
        fieldParams = struct();
    end

    io.assert_atomic_units(sys);

    pos = local_require_field(sys, 'site_pos');
    q   = local_require_field(sys, 'site_charge');
    q   = q(:);

    nSites = size(pos, 1);

    if numel(q) ~= nSites
        error('thole:induced_field_from_charges_periodic:BadChargeSize', ...
            'sys.site_charge must have length N.');
    end

    H = local_get_direct_lattice(sys); %#ok<NASGU>

    exclude_self = true;
    if isfield(fieldParams, 'exclude_self') && ~isempty(fieldParams.exclude_self)
        exclude_self = logical(fieldParams.exclude_self);
    end

    useThole = false;
    if isfield(fieldParams, 'use_thole_damping') && ~isempty(fieldParams.use_thole_damping)
        useThole = logical(fieldParams.use_thole_damping);
    end

    if isfield(fieldParams, 'target_mask') && ~isempty(fieldParams.target_mask)
        target_mask = logical(fieldParams.target_mask(:));
    elseif isfield(sys, 'site_is_polarizable') && ~isempty(sys.site_is_polarizable)
        target_mask = logical(sys.site_is_polarizable(:));
    else
        target_mask = true(nSites, 1);
    end

    if isfield(fieldParams, 'source_mask') && ~isempty(fieldParams.source_mask)
        source_mask = logical(fieldParams.source_mask(:));
    else
        source_mask = abs(q) > 0;
    end

    if numel(target_mask) ~= nSites || numel(source_mask) ~= nSites
        error('thole:induced_field_from_charges_periodic:BadMask', ...
            'target_mask and source_mask must have length N.');
    end

    if ~isfield(fieldParams, 'ewald') || isempty(fieldParams.ewald)
        error('thole:induced_field_from_charges_periodic:MissingEwald', ...
            'fieldParams.ewald is required for periodic mode.');
    end

    ew = fieldParams.ewald;
    local_require_struct_field(ew, 'alpha');
    local_require_struct_field(ew, 'rcut');
    local_require_struct_field(ew, 'kcut');

    boundary = 'tinfoil';
    if isfield(ew, 'boundary') && ~isempty(ew.boundary)
        boundary = lower(string(ew.boundary));
    end
    if boundary ~= "tinfoil"
        error('thole:induced_field_from_charges_periodic:BoundaryNotSupported', ...
            'Only tinfoil boundary is supported in this implementation.');
    end

    if abs(sum(q(source_mask))) > 1e-12
        error('thole:induced_field_from_charges_periodic:NetCharge', ...
            'Periodic source field requires net-neutral selected source charges.');
    end

    src_idx = find(source_mask & (q ~= 0));
    tgt_idx = find(target_mask);

    E = zeros(nSites, 3);
    if isempty(src_idx) || isempty(tgt_idx)
        return;
    end

    % -------------------------
    % real-space cached path
    % -------------------------
    problemRect = struct();
    problemRect.nSites = nSites;
    problemRect.polMask = target_mask;

    scfReal = struct();
    scfReal.use_thole = useThole;

    geomOpts = struct();
    geomOpts.target_mask = target_mask;
    geomOpts.source_mask = source_mask;

    realGeom  = geom.build_periodic_realspace_geom_cache(sys, problemRect, ew, geomOpts);
    realRow   = geom.build_periodic_realspace_row_geom_cache(sys, problemRect, realGeom);
    realCoeff = geom.build_periodic_realspace_charge_coeff_cache(sys, realRow, ew, scfReal);

    Ereal = zeros(numel(tgt_idx), 3);

    row_ptr = realRow.row_ptr;
    srcFull = realRow.source_full_idx;
    dr = realRow.dr;
    coeff_scalar = realCoeff.coeff_scalar(:);

    for a = 1:realRow.nTargetSites
        idx0 = row_ptr(a);
        idx1 = row_ptr(a + 1) - 1;
        if idx1 < idx0
            continue;
        end

        qsrc = q(srcFull(idx0:idx1));
        Ereal(a, :) = -sum((qsrc .* coeff_scalar(idx0:idx1)) .* dr(idx0:idx1, :), 1);
    end

    % -------------------------
    % k-space memory-light path
    % -------------------------
    targetOpts = struct();
    targetOpts.target_mask = target_mask;

    if isfield(fieldParams, 'k_block_size') && ~isempty(fieldParams.k_block_size)
        targetOpts.k_block_size = fieldParams.k_block_size;
    end

    targetCache = geom.build_periodic_kspace_target_cache(sys, problemRect, ew, targetOpts);
    sourceCache = geom.build_periodic_kspace_source_cache(sys, targetCache, source_mask);

    srcCharge = q(sourceCache.sourceSites);
    kCoeff = geom.build_periodic_kspace_charge_coeff_cache(targetCache, sourceCache, srcCharge);

    Erecip = zeros(numel(tgt_idx), 3);

    nk = targetCache.num_kvec;
    if nk > 0
        kBlockSize = targetCache.k_block_size;
        targetPos = targetCache.target_pos;   % nTarget x 3
        kvecs = targetCache.kvecs;            % nk x 3
        pref = kCoeff.pref(:);
        rho_cos = kCoeff.rho_cos(:);
        rho_sin = kCoeff.rho_sin(:);

        for k0 = 1:kBlockSize:nk
            k1 = min(k0 + kBlockSize - 1, nk);
            idx = k0:k1;

            kblk = kvecs(idx, :);                 % nb x 3
            phase = targetPos * kblk.';           % nTarget x nb

            C = cos(phase);
            S = sin(phase);

            % same convention as old code:
            % Sterm = sin_i * rho_cos - cos_i * rho_sin
            Sterm = S .* rho_cos(idx).' - C .* rho_sin(idx).';

            w = Sterm .* pref(idx).';             % nTarget x nb

            Erecip(:, 1) = Erecip(:, 1) + w * kblk(:, 1);
            Erecip(:, 2) = Erecip(:, 2) + w * kblk(:, 2);
            Erecip(:, 3) = Erecip(:, 3) + w * kblk(:, 3);
        end
    end

    E(tgt_idx, :) = Ereal + Erecip;
end

function value = local_require_field(s, name)
    if ~isfield(s, name) || isempty(s.(name))
        error('thole:induced_field_from_charges_periodic:MissingField', ...
            'sys.%s is required.', name);
    end
    value = s.(name);
end

function value = local_require_struct_field(s, name)
    if ~isfield(s, name) || isempty(s.(name))
        error('thole:induced_field_from_charges_periodic:MissingEwaldField', ...
            'fieldParams.ewald.%s is required.', name);
    end
    value = s.(name);
end

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('thole:induced_field_from_charges_periodic:MissingLattice', ...
            'sys.super_lattice or sys.lattice is required.');
    end
end