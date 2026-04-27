function Lmin = shortest_lattice_translation(latOrSysOrH)
%SHORTEST_LATTICE_TRANSLATION Shortest nonzero direct-lattice translation.
%
% Uses project convention:
%   H rows are direct lattice vectors
%   r_cart = f_frac * H

lat = geom.get_lattice(latOrSysOrH);
H = lat.H;

Lmin = inf;
for i = -1:1
    for j = -1:1
        for k = -1:1
            if i == 0 && j == 0 && k == 0
                continue;
            end
            n = [i j k];
            t = n * H;
            Lmin = min(Lmin, norm(t));
        end
    end
end

if ~isfinite(Lmin) || Lmin <= 0
    error('geom:shortest_lattice_translation:Failed', ...
        'Failed to determine a positive shortest lattice translation.');
end
end