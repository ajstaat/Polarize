function X = unstack_xyz(v)
%UNSTACK_XYZ Convert 3N x 1 stacked vector to N x 3 array.

if mod(numel(v), 3) ~= 0
    error('Input length must be a multiple of 3.');
end

X = reshape(v, 3, []).';

end