function idx = block3(i)
%BLOCK3 Return 3-component index block for site i.
%
% Example:
%   idx = util.block3(4)   % returns 10:12

idx = (3*i - 2):(3*i);

end