function v = stack_xyz(X)
%STACK_XYZ Convert N x 3 array to 3N x 1 stacked vector.

if size(X,2) ~= 3
    error('Input must be N x 3.');
end

v = reshape(X.', [], 1);

end