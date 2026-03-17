function Edip = induced_field_from_dipoles(sys, mu, dipoleParams)
%INDUCED_FIELD_FROM_DIPOLES Field at each site from induced dipoles.
%
% Inputs
%   sys           struct with field:
%                 .site_pos   N x 3 Cartesian positions
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
% Formula
%   E(r) = [3 r (mu·r)/r^5] - [mu/r^3]
%
% Optional softening uses:
%   r^2 -> r^2 + softening^2

if nargin < 3 || isempty(dipoleParams)
    dipoleParams = struct();
end

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
    error('sys.site_pos is missing or empty.');
end

pos = sys.site_pos;
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

for a = 1:numel(targetIdx)
    i = targetIdx(a);
    ri = pos(i, :);
    Ei = [0, 0, 0];

    for b = 1:numel(sourceIdx)
        j = sourceIdx(b);

        if exclude_self && i == j
            continue;
        end

        muj = mu(j, :);
        if all(muj == 0)
            continue;
        end

        rij = ri - pos(j, :);
        r2 = dot(rij, rij) + softening^2;

        if r2 == 0
            continue;
        end

        r = sqrt(r2);
        r3 = r2 * r;
        r5 = r2 * r3;

        muDotR = dot(muj, rij);
        Ei = Ei + 3 * rij * (muDotR / r5) - muj / r3;
    end

    Edip(i, :) = Ei;
end

end