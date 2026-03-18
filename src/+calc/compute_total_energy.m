function energy = compute_total_energy(sys, mu, Eext, T)
%COMPUTE_TOTAL_ENERGY Compute induced-dipole energy breakdown.
%
% Inputs
%   sys    working system struct with:
%          .site_is_polarizable
%          .site_alpha
%
%   mu     N x 3 induced dipoles
%   Eext   N x 3 external field
%   T      3N x 3N dipole interaction operator
%
% Output
%   energy struct with fields:
%          .polarization_self
%          .external_charge_dipole
%          .dipole_dipole
%          .total
%
% Energy convention:
%   U_total = U_pol - sum_i mu_i·Eext_i - 1/2 mu^T T mu
%
% where
%   U_pol   = 1/2 sum_i |mu_i|^2 / alpha_i
%   U_ext   = - sum_i mu_i·Eext_i
%   U_dd    = - 1/2 mu^T T mu
%
% Notes
%   - Only polarizable sites contribute to U_pol.
%   - T is assumed to be in the convention:
%         U_pair = 1/2 * mu_vec' * T * mu_vec

if ~isfield(sys, 'site_is_polarizable') || isempty(sys.site_is_polarizable)
    error('sys.site_is_polarizable is missing or empty.');
end
if ~isfield(sys, 'site_alpha') || isempty(sys.site_alpha)
    error('sys.site_alpha is missing or empty.');
end

nSites = sys.n_sites;

if ~isequal(size(mu), [nSites, 3])
    error('mu must be N x 3.');
end
if ~isequal(size(Eext), [nSites, 3])
    error('Eext must be N x 3.');
end
if ~isequal(size(T), [3*nSites, 3*nSites])
    error('T must be 3N x 3N.');
end

polMask = logical(sys.site_is_polarizable(:));
alpha = sys.site_alpha(:);

if numel(alpha) ~= nSites
    error('sys.site_alpha must have length N.');
end

if any(alpha(polMask) <= 0)
    error('Polarizable sites must have positive alpha values.');
end

% 1) Polarization self-cost
mu2 = sum(mu.^2, 2);
Upol = 0.5 * sum(mu2(polMask) ./ alpha(polMask));

% 2) Charge-dipole interaction
Uext = -sum(sum(mu .* Eext, 2));

% 3) Dipole-dipole interaction
mu_vec = util.stack_xyz(mu);
Udd = -0.5 * (mu_vec.' * T * mu_vec);

% Total
Utotal = Upol + Uext + Udd;

energy = struct();
energy.polarization_self = Upol;
energy.external_charge_dipole = Uext;
energy.dipole_dipole = Udd;
energy.total = Utotal;

end