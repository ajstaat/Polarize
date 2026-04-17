function shifts = translation_grid(supercell_size)
%TRANSLATION_GRID Integer translation vectors for a supercell.
%
% Input
%   supercell_size : [nx ny nz]
%
% Output
%   shifts         : (nx*ny*nz) x 3 array of integer shifts
%                    each row is [ix iy iz], with
%                    ix = 0:(nx-1), iy = 0:(ny-1), iz = 0:(nz-1)

    if ~isnumeric(supercell_size) || numel(supercell_size) ~= 3
        error('geom:translation_grid:BadInput', ...
            'supercell_size must be a 1x3 numeric vector.');
    end

    supercell_size = reshape(supercell_size, 1, 3);

    if any(supercell_size < 1) || any(mod(supercell_size,1) ~= 0)
        error('geom:translation_grid:BadSupercellSize', ...
            'supercell_size entries must be positive integers.');
    end

    nx = supercell_size(1);
    ny = supercell_size(2);
    nz = supercell_size(3);

    nTot = nx * ny * nz;
    shifts = zeros(nTot, 3);

    k = 1;
    for ix = 0:(nx-1)
        for iy = 0:(ny-1)
            for iz = 0:(nz-1)
                shifts(k, :) = [ix, iy, iz];
                k = k + 1;
            end
        end
    end
end