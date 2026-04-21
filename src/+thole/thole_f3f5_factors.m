function tf = thole_f3f5_factors(r, alpha_i, alpha_j, a)
%THOLE_F3F5_FACTORS Compute Thole damping factors f3/f5 and corrections l3/l5.
%
% tf = thole.thole_f3f5_factors(r, alpha_i, alpha_j, a)
%
% Inputs
%   r        scalar or array of nonnegative separations
%   alpha_i  scalar or array of nonnegative polarizabilities for site i
%   alpha_j  scalar or array of nonnegative polarizabilities for site j
%   a        scalar or array of nonnegative Thole parameters
%
% Output
%   tf struct with fields:
%     .s
%     .f3
%     .f5
%     .l3
%     .l5
%
% Definitions
%   s  = a * ( r / (alpha_i*alpha_j)^(1/6) )^3
%   f3 = 1 - exp(-s)
%   f5 = 1 - (1 + s) exp(-s)
%   l3 = f3 - 1 = -exp(-s)
%   l5 = f5 - 1 = -(1 + s) exp(-s)
%
% Notes
% - Supports scalar or array inputs.
% - Non-scalar inputs must either be the same size or scalar-expandable.
% - If any of r, alpha_i, alpha_j, or a are zero at an entry, that entry
%   returns the undamped limit:
%       s = 0, f3 = 1, f5 = 1, l3 = 0, l5 = 0

narginchk(4, 4);

validateattributes(r, {'double'}, {'real', 'nonnegative', 'finite'}, ...
    mfilename, 'r', 1);
validateattributes(alpha_i, {'double'}, {'real', 'nonnegative', 'finite'}, ...
    mfilename, 'alpha_i', 2);
validateattributes(alpha_j, {'double'}, {'real', 'nonnegative', 'finite'}, ...
    mfilename, 'alpha_j', 3);
validateattributes(a, {'double'}, {'real', 'nonnegative', 'finite'}, ...
    mfilename, 'a', 4);

[r, alpha_i, alpha_j, a] = local_expand_inputs(r, alpha_i, alpha_j, a);

tf = struct();
tf.s  = zeros(size(r));
tf.f3 = ones(size(r));
tf.f5 = ones(size(r));
tf.l3 = zeros(size(r));
tf.l5 = zeros(size(r));

active = (r > 0) & (alpha_i > 0) & (alpha_j > 0) & (a > 0);
if ~any(active(:))
    return;
end

aij = (alpha_i(active) .* alpha_j(active)).^(1/6);

% Extra guard, though alpha_i/alpha_j > 0 above should already ensure this.
good = (aij > 0);
if ~any(good)
    return;
end

idxActive = find(active);
idxGood = idxActive(good);

s = a(idxGood) .* (r(idxGood) ./ aij(good)).^3;
e = exp(-s);

tf.s(idxGood)  = s;
tf.f3(idxGood) = 1 - e;
tf.f5(idxGood) = 1 - (1 + s) .* e;
tf.l3(idxGood) = -e;
tf.l5(idxGood) = -(1 + s) .* e;

end

function varargout = local_expand_inputs(varargin)
% Expand scalar inputs to match the common non-scalar size.

n = nargin;
sizes = cellfun(@size, varargin, 'UniformOutput', false);
isScalar = cellfun(@isscalar, varargin);

% Find reference size from first non-scalar input.
refSize = [];
for k = 1:n
    if ~isScalar(k)
        refSize = sizes{k};
        break;
    end
end

if isempty(refSize)
    % All scalar: nothing to do.
    varargout = varargin;
    return;
end

% Validate compatibility and expand scalars.
varargout = varargin;
for k = 1:n
    if isScalar(k)
        varargout{k} = repmat(varargin{k}, refSize);
    else
        if ~isequal(size(varargin{k}), refSize)
            error('thole:thole_f3f5_factors:SizeMismatch', ...
                'Inputs must be scalar or all non-scalar inputs must have the same size.');
        end
    end
end
end