function [molA, molB, info] = find_central_molecule_pair(sys, varargin)
%FIND_CENTRAL_MOLECULE_PAIR Choose two molecule images near the supercell center.
%
% Inputs
%   sys working system struct
%
% Name-value options
%   'DifferentBaseMol'   logical, default true
%   'CandidateMolIDs'    numeric vector of allowed unique_mol_id values, default []
%   'CandidateBaseMolIDs' numeric vector of allowed base_mol_id values, default []
%   'TargetSeparation'   numeric scalar, default []
%
% Outputs
%   molA, molB  unique molecule IDs
%   info        struct with center, distances, chosen rows

    p = inputParser;
    addParameter(p, 'DifferentBaseMol', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'CandidateMolIDs', [], @(x) isnumeric(x));
    addParameter(p, 'CandidateBaseMolIDs', [], @(x) isnumeric(x));
    addParameter(p, 'TargetSeparation', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    parse(p, varargin{:});
    opt = p.Results;

    T = builder.list_molecules(sys);

    keep = true(height(T), 1);

    if ~isempty(opt.CandidateMolIDs)
        keep = keep & ismember(T.unique_mol_id, opt.CandidateMolIDs(:));
    end

    if ~isempty(opt.CandidateBaseMolIDs)
        keep = keep & ismember(T.base_mol_id, opt.CandidateBaseMolIDs(:));
    end

    T = T(keep, :);

    if height(T) < 2
        error('Need at least two candidate molecules after filtering.');
    end

    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        center = sum(sys.super_lattice, 1) / 2;
    elseif isfield(sys, 'site_pos') && ~isempty(sys.site_pos)
        mins = min(sys.site_pos, [], 1);
        maxs = max(sys.site_pos, [], 1);
        center = 0.5 * (mins + maxs);
    else
        error('Cannot determine supercell center.');
    end

    xyz = [T.cx, T.cy, T.cz];
    d2center = sqrt(sum((xyz - center).^2, 2));
    T.d2center = d2center;

    [~, orderA] = sort(T.d2center, 'ascend');

    bestScore = inf;
    bestA = [];
    bestB = [];

    for ia = 1:numel(orderA)
        rowA = orderA(ia);

        for ib = 1:numel(orderA)
            rowB = orderA(ib);

            if rowB == rowA
                continue;
            end

            if opt.DifferentBaseMol && T.base_mol_id(rowA) == T.base_mol_id(rowB)
                continue;
            end

            score = T.d2center(rowA) + T.d2center(rowB);

            if ~isempty(opt.TargetSeparation)
                rij = norm(xyz(rowA,:) - xyz(rowB,:));
                score = score + abs(rij - opt.TargetSeparation);
            end

            if score < bestScore
                bestScore = score;
                bestA = rowA;
                bestB = rowB;
            end
        end
    end

    if isempty(bestA) || isempty(bestB)
        error('Could not identify a central molecule pair under the requested constraints.');
    end

    molA = T.unique_mol_id(bestA);
    molB = T.unique_mol_id(bestB);

    info = struct();
    info.center = center;
    info.score = bestScore;
    info.table = T;
    info.rowA = T(bestA, :);
    info.rowB = T(bestB, :);
end