function Edip = induced_field_from_dipoles_thole(sys, mu, dipoleParams)
%INDUCED_FIELD_FROM_DIPOLES_THOLE Field at each site from induced dipoles
%with isotropic Thole damping.
%
% Inputs
%   sys           struct with fields:
%                 .site_pos
%                 .site_alpha
%                 .thole_a
%
%   mu            N x 3 induced dipoles
%
%   dipoleParams  optional struct with fields:
%                 .exclude_self   logical, default true
%                 .softening      scalar, default 0
%                 .target_mask    N x 1 logical, optional
%                 .source_mask    N x 1 logical, optional
%
% Output
%   Edip          N x 3 dipole field at target sites
%
% Bare dipole field tensor:
%   T(r) mu = 3 r (mu·r)/r^5 - mu/r^3
%
% Here we multiply the tensor by a scalar Thole damping factor f(r).

if nargin < 3 || isempty(dipoleParams)
    dipoleParams = struct();
end

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
    error('sys.site_pos is missing or empty.');
end

if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
    error('sys.site_alpha is missing or empty.');
end

if ~isfield(sys, 'thole_a') || isempty(sys.thole_a)
    error('sys.thole_a is missing or empty.');
end

pos = sys.site_pos;
alpha = sys.site_alpha(:);
a = sys.thole_a;

nSites = size(pos, 1);

if ~isequal(size(mu), [nSites, 3])
    error('mu must be N x 3.');
end

exclude_self = true;
if isfield(dipoleParams, 'exclude_self')
    exclude_self = dipoleParams.exclude_self;
end

softening = 0.0;
if isfield(dipoleParams, 'softening')
    softening = dipoleParams.softening;
end

target_mask = true(nSites, 1);
if isfield(dipoleParams, 'target_mask') && ~isempty(dipoleParams.target_mask)
    target_mask = logical(dipoleParams.target_mask(:));
end

source_mask = true(nSites, 1);
if isfield(dipoleParams, 'source_mask') && ~isempty(dipoleParams.source_mask)
    source_mask = logical(dipoleParams.source_mask(:));
end

Edip = zeros(nSites, 3);

targetIdx = find(target_mask);
sourceIdx = find(source_mask);

for aa = 1:numel(targetIdx)
    i = targetIdx(aa);
    ri = pos(i, :);
    Ei = [0, 0, 0];

    for bb = 1:numel(sourceIdx)
        j = sourceIdx(bb);

        if exclude_self && i == j
            continue;
        end

        muj = mu(j, :);
        if all(muj == 0)
            continue;
        end

        rij = ri - pos(j, :);
        r2_bare = dot(rij, rij);

        if r2_bare == 0
            continue;
        end

        r2 = r2_bare + softening^2;
        r = sqrt(r2);
        r3 = r2 * r;
        r5 = r2 * r3;

        muDotR = dot(muj, rij);

        % Bare dipole field
        Ebare = 3 * rij * (muDotR / r5) - muj / r3;

        % Thole damping factor based on physical separation
        r_phys = sqrt(r2_bare);
        f = thole.thole_scale_factor(r_phys, alpha(i), alpha(j), a);

        Ei = Ei + f * Ebare;
    end

    Edip(i, :) = Ei;
end

end