function [Urecip, meta] = reciprocal_space_energy_triclinic(mu, r, H, alpha, kcut)
%RECIPROCAL_SPACE_ENERGY_TRICLINIC Reciprocal-space part of triclinic dipole Ewald energy.
%
% [Urecip, meta] = ewald.reciprocal_space_energy_triclinic(mu, r, H, alpha, kcut)
%
% Inputs
%   mu     N x 3 dipole vectors
%   r      N x 3 Cartesian positions
%   H      3 x 3 direct lattice matrix, columns are lattice vectors
%   alpha  Ewald screening parameter
%   kcut   reciprocal-space cutoff in |k|
%
% Outputs
%   Urecip reciprocal-space contribution to dipole Ewald energy
%   meta   struct with bookkeeping info
%          .num_kvec
%          .k2
%          .knorm
%
% Notes
%   - Electrostatic prefactor is 1.
%   - Uses the structure-factor form:
%       Urecip = sum_k pref(k) * |sum_j (k·mu_j) exp(i k·r_j)|^2

    validateattributes(mu, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
        mfilename, 'mu', 1);
    validateattributes(r, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
        mfilename, 'r', 2);
    if size(mu, 1) ~= size(r, 1)
        error('ewald:reciprocal_space_energy_triclinic:SizeMismatch', ...
            'mu and r must have the same number of rows.');
    end
    validateattributes(H, {'double'}, {'size', [3, 3], 'finite', 'real'}, ...
        mfilename, 'H', 3);
    validateattributes(alpha, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'alpha', 4);
    validateattributes(kcut, {'double'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'kcut', 5);

    V = abs(det(H));
    if V <= 1e-14
        error('ewald:reciprocal_space_energy_triclinic:SingularCell', ...
            'Cell matrix must have positive nonzero volume.');
    end

    [kvecs, kmeta] = ewald.enumerate_kvecs_triclinic(H, kcut);

    nk = size(kvecs, 1);
    if nk == 0
        Urecip = 0.0;
        meta = struct();
        meta.num_kvec = 0;
        meta.k2 = zeros(0,1);
        meta.knorm = zeros(0,1);
        return;
    end

    k2 = kmeta.k2;
    pref = (2 * pi / V) * exp(-k2 / (4 * alpha^2)) ./ k2;

    Urecip = 0.0;

    for m = 1:nk
        kvec = kvecs(m, :);

        kdotmu = mu * kvec.';      % N x 1
        phase = r * kvec.';        % N x 1
        Smu = sum(kdotmu .* exp(1i * phase));

        Urecip = Urecip + pref(m) * abs(Smu)^2;
    end

    meta = struct();
    meta.num_kvec = nk;
    meta.k2 = k2;
    meta.knorm = kmeta.knorm;
end