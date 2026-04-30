function E = apply_periodic_real_cache(realCache, mu)
%APPLY_PERIODIC_REAL_CACHE Apply cached periodic real-space dipole operator.
%
% E = thole.apply_periodic_real_cache(realCache, mu)
%
% Fast path: assumes realCache has already been normalized by the operator
% builder. Do expensive validation/conversion at build time, not at every
% GMRES/Jacobi matvec.

nPol = size(mu, 1);
E = zeros(nPol, 3);

if realCache.n_entries == 0
    return;
end

if isfield(realCache, 'row_ind') && ~isempty(realCache.row_ind)
    rowInd = realCache.row_ind;
else
    rowInd = local_row_indices_from_row_ptr(realCache.row_ptr);
end

cols = realCache.col_idx;
dr = realCache.dr;

muSrc = mu(cols, :);
muDotR = sum(muSrc .* dr, 2);

Eentry = realCache.coeff_iso .* muSrc + ...
    realCache.coeff_dyad .* muDotR .* dr;

E(:, 1) = accumarray(rowInd, Eentry(:, 1), [nPol 1], @sum, 0);
E(:, 2) = accumarray(rowInd, Eentry(:, 2), [nPol 1], @sum, 0);
E(:, 3) = accumarray(rowInd, Eentry(:, 3), [nPol 1], @sum, 0);
end

function rowInd = local_row_indices_from_row_ptr(rowPtr)
counts = diff(rowPtr(:));
nRows = numel(counts);
rowInd = repelem((1:nRows).', counts);
end