function sCache = build_periodic_kspace_source_cache(sys, targetCache, sourceMask)
%BUILD_PERIODIC_KSPACE_SOURCE_CACHE Build source-side phase cache reusing target k-geometry.
%
% sCache = geom.build_periodic_kspace_source_cache(sys, targetCache, sourceMask)
%
% Inputs
%   sys         canonical polarization-system struct
%   targetCache struct from geom.build_periodic_kspace_target_cache(...)
%   sourceMask  N x 1 logical
%
% Output
%   sCache struct with fields:
%       .mode
%       .nSites
%       .source_mask
%       .sourceSites
%       .same_as_target
%       .source_cos
%       .source_sin

    narginchk(3, 3);

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('geom:build_periodic_kspace_source_cache:MissingSitePos', ...
            'sys.site_pos is required.');
    end

    pos = sys.site_pos;
    nSites = size(pos, 1);

    sourceMask = logical(sourceMask(:));
    if numel(sourceMask) ~= nSites
        error('geom:build_periodic_kspace_source_cache:BadSourceMask', ...
            'sourceMask must have length nSites.');
    end

    sourceSites = find(sourceMask);
    sameAsTarget = isequal(sourceMask, targetCache.target_mask);

    if sameAsTarget
        source_cos = targetCache.target_cos;
        source_sin = targetCache.target_sin;
    else
        nk = targetCache.num_kvec;
        if nk == 0
            source_cos = zeros(numel(sourceSites), 0);
            source_sin = zeros(numel(sourceSites), 0);
        else
            pos_source = pos(sourceSites, :);
            phase_source = pos_source * targetCache.kvecs.';
            source_cos = cos(phase_source);
            source_sin = sin(phase_source);
        end
    end

    sCache = struct();
    sCache.mode = 'periodic_kspace_source';
    sCache.nSites = nSites;
    sCache.source_mask = sourceMask;
    sCache.sourceSites = sourceSites;
    sCache.same_as_target = sameAsTarget;
    sCache.source_cos = source_cos;
    sCache.source_sin = source_sin;
end