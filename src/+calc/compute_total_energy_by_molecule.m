function Emol = compute_total_energy_by_molecule(sys, mu, Eext)
%COMPUTE_TOTAL_ENERGY_BY_MOLECULE Per-molecule local energy contributions.
%
% Inputs
%   sys    with fields:
%          .site_mol_id
%          .site_is_polarizable
%          .site_alpha
%
%   mu     N x 3
%   Eext   N x 3
%
% Output
%   Emol   table with columns:
%          unique_mol_id
%          polarization_self
%          external_charge_dipole
%          subtotal_local
%
% Notes
%   - This routine only includes local terms:
%         U_pol and U_ext
%   - It does not yet partition the pairwise dipole-dipole term.

if ~isfield(sys, 'site_mol_id') || isempty(sys.site_mol_id)
    error('sys.site_mol_id is missing or empty.');
end
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

molIDs = unique(sys.site_mol_id(:));
polMask = logical(sys.site_is_polarizable(:));
alpha = sys.site_alpha(:);

nMol = numel(molIDs);
Upol = zeros(nMol, 1);
Uext = zeros(nMol, 1);

for m = 1:nMol
    molID = molIDs(m);
    idx = find(sys.site_mol_id == molID);

    idxPol = idx(polMask(idx));
    if ~isempty(idxPol)
        mu2 = sum(mu(idxPol, :).^2, 2);
        Upol(m) = 0.5 * sum(mu2 ./ alpha(idxPol));
    end

    Uext(m) = -sum(sum(mu(idx, :) .* Eext(idx, :), 2));
end

Emol = table(molIDs, Upol, Uext, Upol + Uext, ...
    'VariableNames', {'unique_mol_id', 'polarization_self', 'external_charge_dipole', 'subtotal_local'});
end