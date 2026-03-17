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

if numel(supercell_size) ~= 3
    error('supercell_size must be a 1x3 vector.');
end

nx = supercell_size(1);
ny = supercell_size(2);
nz = supercell_size(3);

if any(supercell_size < 1) || any(mod(supercell_size,1) ~= 0)
    error('supercell_size entries must be positive integers.');
end

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