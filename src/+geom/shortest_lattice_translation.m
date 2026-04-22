function Lmin = shortest_lattice_translation(H)
%SHORTEST_LATTICE_TRANSLATION Shortest nonzero direct-lattice translation.
%
%   Lmin = geom.shortest_lattice_translation(H)
%
% Input
%   H : 3x3 direct lattice matrix with lattice vectors as COLUMNS
%
% Output
%   Lmin : shortest nonzero lattice translation length
%
% Notes
%   - For the single-image condition used in periodic real-space work,
%     the relevant quantity is the norm of the shortest nonzero lattice
%     vector H*n for integer n ~= 0.
%   - In practice, checking n in {-1,0,1}^3 \ {0} is sufficient for the
%     shortest translation in ordinary unit-cell geometry.

validateattributes(H, {'double'}, {'size',[3 3], 'finite', 'real'}, ...
    mfilename, 'H', 1);

Lmin = inf;

for i = -1:1
    for j = -1:1
        for k = -1:1
            if i == 0 && j == 0 && k == 0
                continue;
            end

            t = H * [i; j; k];
            Lmin = min(Lmin, norm(t));
        end
    end
end

if ~isfinite(Lmin) || Lmin <= 0
    error('geom:shortest_lattice_translation:Failed', ...
        'Failed to determine a positive shortest lattice translation.');
end
end