function Tij = reciprocal_space_tensor_block_triclinic(ri, rj, H, alpha, kcut, kvecs)
%RECIPROCAL_SPACE_TENSOR_BLOCK_TRICLINIC Reciprocal-space Ewald dipole tensor block.
%
% Inputs
%   ri, rj   1x3 Cartesian positions
%   H        3x3 direct lattice matrix, columns are lattice vectors
%   alpha    Ewald screening parameter
%   kcut     reciprocal-space cutoff
%   kvecs    optional Nk x 3 array from ewald.enumerate_kvecs_triclinic
%
% Output
%   Tij      3x3 reciprocal-space tensor block such that
%            mu_i' * Tij * mu_j contributes to reciprocal energy
%
% Formula
%   Tij = sum_{k != 0} (2*pi/V) exp(-k^2/(4 alpha^2)) / k^2
%         * cos(k·(ri-rj)) * (k k^T)

if numel(ri) ~= 3 || numel(rj) ~= 3
    error('ri and rj must be 1x3 vectors.');
end
if ~isequal(size(H), [3 3])
    error('H must be 3x3.');
end
if alpha <= 0 || kcut <= 0
    error('alpha and kcut must be positive.');
end

V = abs(det(H));
if V <= 0
    error('Cell matrix must have positive nonzero volume.');
end

if nargin < 6 || isempty(kvecs)
    [kvecs, ~] = ewald.enumerate_kvecs_triclinic(H, kcut);
end

rij = ri - rj;
Tij = zeros(3,3);

for m = 1:size(kvecs,1)
    kvec = kvecs(m,:);
    k2 = dot(kvec, kvec);

    pref = (2*pi / V) * exp(-k2 / (4*alpha^2)) / k2;
    phase = dot(kvec, rij);

    kkT = kvec(:) * kvec(:).';
    Tij = Tij + 2 * pref * cos(phase) * kkT;
end

end