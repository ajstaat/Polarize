function E = induced_field_from_charges(sys, fieldParams)
%INDUCED_FIELD_FROM_CHARGES Electric field at each site from point charges.
%
% Inputs
%   sys          working system struct with fields:
%                .site_pos     N x 3 Cartesian positions
%                .site_charge  N x 1 charges
%
%   fieldParams  struct, optional fields:
%                .exclude_self           logical, default true
%                .softening              scalar, default 0
%                .target_mask            N x 1 logical, optional
%                .source_mask            N x 1 logical, optional
%
% Output
%   E            N x 3 electric field at each target site
%
% Notes
%   Field formula used:
%       E_i = sum_j q_j * r_ij / |r_ij|^3
%   with optional softening:
%       |r_ij|^2 -> |r_ij|^2 + softening^2

if nargin < 2 || isempty(fieldParams)
    fieldParams = struct();
end

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
    error('sys.site_pos is missing or empty.');
end

if ~isfield(sys, 'site_charge') || isempty(sys.site_charge)
    error('sys.site_charge is missing or empty.');
end

pos = sys.site_pos;
q = sys.site_charge(:);

nSites = size(pos, 1);

exclude_self = true;
if isfield(fieldParams, 'exclude_self')
    exclude_self = fieldParams.exclude_self;
end

softening = 0.0;
if isfield(fieldParams, 'softening')
    softening = fieldParams.softening;
end

target_mask = true(nSites, 1);
if isfield(fieldParams, 'target_mask') && ~isempty(fieldParams.target_mask)
    target_mask = logical(fieldParams.target_mask(:));
end

source_mask = true(nSites, 1);
if isfield(fieldParams, 'source_mask') && ~isempty(fieldParams.source_mask)
    source_mask = logical(fieldParams.source_mask(:));
end

if numel(target_mask) ~= nSites
    error('fieldParams.target_mask must have length N.');
end

if numel(source_mask) ~= nSites
    error('fieldParams.source_mask must have length N.');
end

E = zeros(nSites, 3);

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

        qj = q(j);
        if qj == 0
            continue;
        end

        rij = ri - pos(j, :);
        r2 = dot(rij, rij) + softening^2;

        if r2 == 0
            continue;
        end

        r3 = r2^(3/2);
        Ei = Ei + qj * rij / r3;
    end

    E(i, :) = Ei;
end

end