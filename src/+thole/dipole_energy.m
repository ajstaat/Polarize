function Epol = dipole_energy(sys, mu, Eext)
%DIPOLE_ENERGY Simple polarization energy diagnostic.
%
% Inputs
%   sys    with .site_is_polarizable
%   mu     N x 3 induced dipoles
%   Eext   N x 3 external field
%
% Output
%   Epol   scalar
%
% Uses the common expression:
%   Epol = -1/2 sum_i mu_i · Eext_i
%
% This is a useful first-stage diagnostic. Later, once the full periodic
% interaction model is assembled, you may replace this with a more complete
% energy expression.

if ~isfield(sys, 'site_is_polarizable') || isempty(sys.site_is_polarizable)
    error('sys.site_is_polarizable is missing or empty.');
end

polMask = logical(sys.site_is_polarizable(:));

dotVals = sum(mu(polMask, :) .* Eext(polMask, :), 2);
Epol = -0.5 * sum(dotVals);

end