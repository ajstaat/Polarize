function [Ereal, info] = realspace_charge_field_ewald(sys, opts)
%REALSPACE_CHARGE_FIELD_EWALD Direct real-space Ewald charge field.
%
% Computes the short-range real-space contribution to the periodic Ewald
% charge field:
%
%   E_i = sum_{j,R} q_j r_ijR *
%         [ erfc(alpha r)/r^3
%           + 2 alpha/sqrt(pi) exp(-alpha^2 r^2)/r^2 ]
%
% Inputs
%   sys.site_pos
%   sys.site_charge
%   sys.site_is_polarizable
%   sys.super_lattice or sys.lattice
%
% opts.alpha
% opts.rcut
% opts.target_mask
% opts.source_mask
% opts.exclude_self default true
%
% Output
%   Ereal  N x 3 field, non-target rows zero

    io.assert_atomic_units(sys);

    alpha = opts.alpha;
    rcut = opts.rcut;

    pos = sys.site_pos;
    q = sys.site_charge(:);
    nSites = size(pos, 1);

    targetMask = local_get_mask(opts, 'target_mask', true(nSites,1));
    sourceMask = local_get_mask(opts, 'source_mask', abs(q) > 0);

    excludeSelf = local_get_opt(opts, 'exclude_self', true);

    targetSites = find(targetMask);
    sourceSites = find(sourceMask);

    Hrow = local_get_direct_lattice(sys);
    Hcol = Hrow.';

    Ereal = zeros(nSites, 3);

    vecLens = vecnorm(Hcol, 2, 1);
    nMax = ceil(rcut / min(vecLens)) + 2;

    nTerms = 0;

    for ia = 1:numel(targetSites)
        i = targetSites(ia);
        ri = pos(i, :);

        Ei = [0.0, 0.0, 0.0];

        for n1 = -nMax:nMax
            for n2 = -nMax:nMax
                for n3 = -nMax:nMax
                    R = (Hcol * [n1; n2; n3]).';

                    rj = pos(sourceSites, :) + R;
                    dr = ri - rj;
                    r2 = sum(dr.^2, 2);

                    keep = r2 > 1e-24 & r2 <= rcut^2;

                    if excludeSelf && n1 == 0 && n2 == 0 && n3 == 0
                        keep(sourceSites == i) = false;
                    end

                    if ~any(keep)
                        continue;
                    end

                    srcKeep = sourceSites(keep);
                    drk = dr(keep, :);
                    r2k = r2(keep);
                    rk = sqrt(r2k);

                    pref = erfc(alpha .* rk) ./ (rk.^3) + ...
                        (2 * alpha / sqrt(pi)) .* exp(-(alpha^2) .* r2k) ./ r2k;

                    Ei = Ei + sum((q(srcKeep) .* pref) .* drk, 1);

                    nTerms = nTerms + nnz(keep);
                end
            end
        end

        Ereal(i, :) = Ei;
    end

    info = struct();
    info.nTargets = numel(targetSites);
    info.nSources = numel(sourceSites);
    info.nTerms = nTerms;
    info.alpha = alpha;
    info.rcut = rcut;
end

function H = local_get_direct_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('p3m:realspace_charge_field_ewald:MissingLattice', ...
            'Missing sys.super_lattice or sys.lattice.');
    end
end

function mask = local_get_mask(opts, name, defaultMask)
    if isfield(opts, name) && ~isempty(opts.(name))
        mask = logical(opts.(name)(:));
    else
        mask = logical(defaultMask(:));
    end
end

function val = local_get_opt(s, name, defaultVal)
    if isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = defaultVal;
    end
end