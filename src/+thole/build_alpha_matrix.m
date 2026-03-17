function A = build_alpha_matrix(sys)
%BUILD_ALPHA_MATRIX Build 3N x 3N block-diagonal polarizability matrix.
%
% Input
%   sys    struct with fields:
%          .site_alpha
%          .site_is_polarizable
%
% Output
%   A      3N x 3N block-diagonal matrix
%
% For isotropic scalar polarizabilities:
%   block i = alpha_i * I3 for polarizable sites
%   block i = 0 otherwise

if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
    error('sys.site_alpha is missing or empty.');
end

if ~isfield(sys, 'site_is_polarizable') || isempty(sys.site_is_polarizable)
    error('sys.site_is_polarizable is missing or empty.');
end

nSites = sys.n_sites;
alpha = sys.site_alpha(:);
polMask = logical(sys.site_is_polarizable(:));

if numel(alpha) ~= nSites
    error('sys.site_alpha must have length N.');
end

A = zeros(3*nSites, 3*nSites);
I3 = eye(3);

for i = 1:nSites
    if ~polMask(i)
        continue;
    end

    ii = util.block3(i);
    A(ii, ii) = alpha(i) * I3;
end

end