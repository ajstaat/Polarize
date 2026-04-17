function [refMolID, summary] = choose_center_reference_molecule(sys, varargin)
%CHOOSE_CENTER_REFERENCE_MOLECULE Choose the complete molecule nearest the
% geometric center of the displayed supercell.
%
% [refMolID, summary] = builder.choose_center_reference_molecule(sys)
% [refMolID, summary] = builder.choose_center_reference_molecule(sys, ...)
%
% Optional name-value inputs
%   'RequireComplete'  logical, default true
%   'Verbose'          logical, default false
%
% Output
%   refMolID   scalar supercell molecule ID
%
%   summary    struct with fields:
%       .supercell_center
%       .candidate_mol_ids
%       .candidate_com
%       .candidate_distance
%       .chosen_mol_id
%       .chosen_com

    p = inputParser;
    addRequired(p, 'sys', @isstruct);
    addParameter(p, 'RequireComplete', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
    parse(p, sys, varargin{:});
    opt = p.Results;

    validate_sys(sys);

    T = sys.molecule_table;
    center = 0.5 * sum(sys.super_lattice, 1);

    nTotal = numel(T.molecule_id);
    nComplete = nnz(T.is_complete_in_display);

    if opt.RequireComplete
        rows = find(T.is_complete_in_display);
        if isempty(rows)
            error('builder:choose_center_reference_molecule:NoCompleteMolecules', ...
                'No complete molecules are available in the displayed supercell.');
        end
    else
        rows = (1:numel(T.molecule_id)).';
    end

    molIDs = T.molecule_id(rows);
    com = T.com(rows, :);

    d2 = sum((com - center).^2, 2);
    [~, iMin] = min(d2);

    refMolID = molIDs(iMin);

    summary = struct();
    summary.supercell_center = center;
    summary.candidate_mol_ids = molIDs;
    summary.candidate_com = com;
    summary.candidate_distance = sqrt(d2);
    summary.chosen_mol_id = refMolID;
    summary.chosen_com = com(iMin, :);

    if opt.Verbose
        fprintf('Center-reference selection summary:\n');
        fprintf('  total molecules found        = %d\n', nTotal);
        fprintf('  complete molecules found     = %d\n', nComplete);
        fprintf('  candidates considered        = %d\n', numel(molIDs));
        fprintf('  chosen reference ID          = %d\n', refMolID);
        fprintf('  chosen COM                   = [%9.4f %9.4f %9.4f]\n', summary.chosen_com);
        fprintf('  supercell center             = [%9.4f %9.4f %9.4f]\n', center);
        fprintf('  distance to center           = %9.4f\n', summary.candidate_distance(iMin));

        if numel(molIDs) <= 12
            fprintf('  candidate IDs                = %s\n', mat2str(molIDs(:).'));
        end
    end
end


function validate_sys(sys)
    if ~isstruct(sys)
        error('builder:choose_center_reference_molecule:BadInput', ...
            'sys must be a struct.');
    end

    if ~isfield(sys, 'super_lattice') || isempty(sys.super_lattice)
        error('builder:choose_center_reference_molecule:MissingSuperLattice', ...
            'sys.super_lattice is required and missing/empty.');
    end

    if ~isequal(size(sys.super_lattice), [3 3])
        error('builder:choose_center_reference_molecule:BadSuperLattice', ...
            'sys.super_lattice must be 3 x 3.');
    end

    if ~isfield(sys, 'molecule_table') || isempty(sys.molecule_table)
        error('builder:choose_center_reference_molecule:MissingMoleculeTable', ...
            'sys.molecule_table is required and missing/empty.');
    end

    T = sys.molecule_table;
    required = {'molecule_id', 'com', 'is_complete_in_display'};

    for k = 1:numel(required)
        name = required{k};
        if ~isfield(T, name) || isempty(T.(name))
            error('builder:choose_center_reference_molecule:BadMoleculeTable', ...
                'sys.molecule_table.%s is required and missing/empty.', name);
        end
    end

    if size(T.com, 2) ~= 3
        error('builder:choose_center_reference_molecule:BadCOM', ...
            'sys.molecule_table.com must be N x 3.');
    end
end