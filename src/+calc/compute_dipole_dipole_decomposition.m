function dd = compute_dipole_dipole_decomposition(sys, mu, T)
%COMPUTE_DIPOLE_DIPOLE_DECOMPOSITION Decompose dipole-dipole energy by active/environment groups.
%
% Inputs
%   sys   working system struct with fields:
%         .site_mol_id
%         .active_molecules
%
%   mu    N x 3 induced dipoles
%   T     3N x 3N dipole interaction operator
%
% Output
%   dd    struct with fields:
%         .active_active
%         .active_environment
%         .environment_environment
%         .total
%
% Energy convention:
%   U_dd = -1/2 * mu_vec' * T * mu_vec
%
% Group decomposition:
%   U_dd(AA) = -1/2 * mu_A' T_AA mu_A
%   U_dd(AE) = -      mu_A' T_AE mu_E
%   U_dd(EE) = -1/2 * mu_E' T_EE mu_E

if ~isfield(sys, 'site_mol_id') || isempty(sys.site_mol_id)
    error('sys.site_mol_id is missing or empty.');
end
if ~isfield(sys, 'active_molecules')
    error('sys.active_molecules field is missing.');
end

nSites = sys.n_sites;

if ~isequal(size(mu), [nSites, 3])
    error('mu must be N x 3.');
end
if ~isequal(size(T), [3*nSites, 3*nSites])
    error('T must be 3N x 3N.');
end

activeMolIDs = sys.active_molecules(:).';
siteActiveMask = ismember(sys.site_mol_id(:), activeMolIDs);
siteEnvMask = ~siteActiveMask;

idxA = [];
idxE = [];

for i = 1:nSites
    ii = util.block3(i);
    if siteActiveMask(i)
        idxA = [idxA, ii]; %#ok<AGROW>
    else
        idxE = [idxE, ii]; %#ok<AGROW>
    end
end

mu_vec = util.stack_xyz(mu);

if isempty(idxA)
    Uaa = 0.0;
    Uae = 0.0;
else
    muA = mu_vec(idxA);
    TAA = T(idxA, idxA);
    Uaa = -0.5 * (muA.' * TAA * muA);

    if isempty(idxE)
        Uae = 0.0;
    else
        muE = mu_vec(idxE);
        TAE = T(idxA, idxE);
        Uae = -(muA.' * TAE * muE);
    end
end

if isempty(idxE)
    Uee = 0.0;
else
    muE = mu_vec(idxE);
    TEE = T(idxE, idxE);
    Uee = -0.5 * (muE.' * TEE * muE);
end

dd = struct();
dd.active_active = Uaa;
dd.active_environment = Uae;
dd.environment_environment = Uee;
dd.total = Uaa + Uae + Uee;

end