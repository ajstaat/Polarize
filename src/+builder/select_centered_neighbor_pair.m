function result = select_centered_neighbor_pair(sys, varargin)
%SELECT_CENTERED_NEIGHBOR_PAIR Select the best-centered pair for a requested
% neighbor relation by searching over complete reference molecules.
%
% result = builder.select_centered_neighbor_pair(sys, ...)
%
% Required/important name-value inputs
%   'Relation'     : 'same_stack' | 'side_stack' | 'skew'
%
% Relation-specific required inputs
%   For 'same_stack':
%       'Shell'     : positive integer
%
%   For 'side_stack':
%       'Shell'     : positive integer
%       'Member'    : positive integer
%
%   For 'skew':
%       'Shell'     : positive integer
%
% Optional name-value inputs
%   'StackAxis'         : 'a' | 'b' | 'c' | numeric 1x3 vector, default 'b'
%   'Direction'         : 'either' | '+' | '-', default 'either'
%
%   Same-stack options
%   'PerpTol'           : default 1e-3
%   'ShellTol'          : default 1e-3
%
%   Side-stack options
%   'SameStackPerpTol'  : default 1e-3
%   'VectorGroupTol'    : default 1e-3
%
%   Skew options
%   'NormalTargetDeg'   : default 72
%   'NormalTolDeg'      : default 15
%
%   Common
%   'Verbose'           : logical, default false
%
% Output
%   result struct with fields:
%       .relation
%       .reference_mol_id
%       .neighbor_mol_id
%       .reference_com
%       .neighbor_com
%       .pair_midpoint
%       .supercell_center
%       .midpoint_distance
%       .reference_desc
%       .selector_result
%       .candidate_table
%
% candidate_table columns
%   reference_mol_id
%   neighbor_mol_id
%   ref_com_x, ref_com_y, ref_com_z
%   nbr_com_x, nbr_com_y, nbr_com_z
%   mid_x, mid_y, mid_z
%   midpoint_distance

    p = inputParser;
    addRequired(p, 'sys', @isstruct);

    addParameter(p, 'Relation', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Shell', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1 && mod(x,1)==0));
    addParameter(p, 'Member', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1 && mod(x,1)==0));

    addParameter(p, 'StackAxis', 'b', @(x) ischar(x) || isstring(x) || (isnumeric(x) && numel(x) == 3));
    addParameter(p, 'Direction', 'either', @(x) ischar(x) || isstring(x));

    addParameter(p, 'PerpTol', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'ShellTol', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'SameStackPerpTol', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'VectorGroupTol', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'NormalTargetDeg', 72, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'NormalTolDeg', 15, @(x) isnumeric(x) && isscalar(x) && x >= 0);

    addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
    parse(p, sys, varargin{:});
    opt = p.Results;

    relation = char(string(opt.Relation));
    if ~ismember(relation, {'same_stack','side_stack','skew'})
        error('builder:select_centered_neighbor_pair:BadRelation', ...
            'Relation must be ''same_stack'', ''side_stack'', or ''skew''.');
    end

    if isempty(opt.Shell)
        error('builder:select_centered_neighbor_pair:MissingShell', ...
            'Shell is required.');
    end

    if strcmp(relation, 'side_stack') && isempty(opt.Member)
        error('builder:select_centered_neighbor_pair:MissingMember', ...
            'Member is required for Relation = ''side_stack''.');
    end

    validate_sys(sys);

    completeIDs = builder.complete_molecule_ids(sys);
    if isempty(completeIDs)
        error('builder:select_centered_neighbor_pair:NoCompleteMolecules', ...
            'No complete molecules are available in the displayed supercell.');
    end

    center = 0.5 * sum(sys.super_lattice, 1);

    refIDs = [];
    nbrIDs = [];
    refCOMs = zeros(0,3);
    nbrCOMs = zeros(0,3);
    mids = zeros(0,3);
    midDist = zeros(0,1);
    descList = {};
    selList = {};

    for iRef = 1:numel(completeIDs)
        refMolID = completeIDs(iRef);

        desc = builder.complete_molecule_descriptors_relative_to_reference(sys, refMolID, ...
            'StackAxis', opt.StackAxis, ...
            'IncludeNormals', true, ...
            'Verbose', false);

        try
            switch relation
                case 'same_stack'
                    sel = builder.select_same_stack_neighbor(desc, opt.Shell, ...
                        'PerpTol', opt.PerpTol, ...
                        'Direction', opt.Direction, ...
                        'ShellTol', opt.ShellTol, ...
                        'Verbose', false);

                case 'side_stack'
                    sel = builder.select_side_stack_neighbor(desc, opt.Shell, opt.Member, ...
                        'SameStackPerpTol', opt.SameStackPerpTol, ...
                        'Direction', opt.Direction, ...
                        'ShellTol', opt.ShellTol, ...
                        'VectorGroupTol', opt.VectorGroupTol, ...
                        'Verbose', false);

                case 'skew'
                    sel = builder.select_skew_neighbor(desc, opt.Shell, ...
                        'NormalTargetDeg', opt.NormalTargetDeg, ...
                        'NormalTolDeg', opt.NormalTolDeg, ...
                        'Direction', opt.Direction, ...
                        'ShellTol', opt.ShellTol, ...
                        'Verbose', false);
            end
        catch
            continue;
        end

        if isempty(sel.match_table)
            continue;
        end

        refIdx = find(sys.molecule_table.molecule_id == refMolID, 1, 'first');
        thisRefCOM = sys.molecule_table.com(refIdx, :);

        for j = 1:height(sel.match_table)
            nbrMolID = sel.match_table.molecule_id(j);

            % Avoid double counting exact reversed duplicate pairs later
            nbrIdx = find(sys.molecule_table.molecule_id == nbrMolID, 1, 'first');
            thisNbrCOM = sys.molecule_table.com(nbrIdx, :);

            mid = 0.5 * (thisRefCOM + thisNbrCOM);
            dmid = norm(mid - center);

            refIDs(end+1,1) = refMolID; %#ok<AGROW>
            nbrIDs(end+1,1) = nbrMolID; %#ok<AGROW>
            refCOMs(end+1,:) = thisRefCOM; %#ok<AGROW>
            nbrCOMs(end+1,:) = thisNbrCOM; %#ok<AGROW>
            mids(end+1,:) = mid; %#ok<AGROW>
            midDist(end+1,1) = dmid; %#ok<AGROW>
            descList{end+1,1} = desc; %#ok<AGROW>
            selList{end+1,1} = sel; %#ok<AGROW>
        end
    end

    if isempty(refIDs)
        error('builder:select_centered_neighbor_pair:NoValidPairs', ...
            'No valid complete pairs were found for the requested relation.');
    end

    candidateTable = table( ...
        refIDs, nbrIDs, ...
        refCOMs(:,1), refCOMs(:,2), refCOMs(:,3), ...
        nbrCOMs(:,1), nbrCOMs(:,2), nbrCOMs(:,3), ...
        mids(:,1), mids(:,2), mids(:,3), ...
        midDist, ...
        'VariableNames', { ...
            'reference_mol_id', 'neighbor_mol_id', ...
            'ref_com_x', 'ref_com_y', 'ref_com_z', ...
            'nbr_com_x', 'nbr_com_y', 'nbr_com_z', ...
            'mid_x', 'mid_y', 'mid_z', ...
            'midpoint_distance'});

    [~, order] = sort(candidateTable.midpoint_distance, 'ascend');
    candidateTable = candidateTable(order, :);
    descList = descList(order);
    selList = selList(order);

    result = struct();
    result.relation = relation;
    result.reference_mol_id = candidateTable.reference_mol_id(1);
    result.neighbor_mol_id = candidateTable.neighbor_mol_id(1);
    result.reference_com = [candidateTable.ref_com_x(1), candidateTable.ref_com_y(1), candidateTable.ref_com_z(1)];
    result.neighbor_com = [candidateTable.nbr_com_x(1), candidateTable.nbr_com_y(1), candidateTable.nbr_com_z(1)];
    result.pair_midpoint = [candidateTable.mid_x(1), candidateTable.mid_y(1), candidateTable.mid_z(1)];
    result.supercell_center = center;
    result.midpoint_distance = candidateTable.midpoint_distance(1);
    result.reference_desc = descList{1};
    result.selector_result = selList{1};
    result.candidate_table = candidateTable;

    if opt.Verbose
        fprintf('Centered-pair selector summary:\n');
        fprintf('  relation                = %s\n', relation);
        fprintf('  complete references     = %d\n', numel(completeIDs));
        fprintf('  valid pairs found       = %d\n', height(candidateTable));
        fprintf('  chosen reference ID     = %d\n', result.reference_mol_id);
        fprintf('  chosen neighbor ID      = %d\n', result.neighbor_mol_id);
        fprintf('  midpoint distance       = %.4f\n', result.midpoint_distance);
        fprintf('  pair midpoint           = [%9.4f %9.4f %9.4f]\n', result.pair_midpoint);
        fprintf('  supercell center        = [%9.4f %9.4f %9.4f]\n', center);
    end
end


function validate_sys(sys)
    if ~isstruct(sys)
        error('builder:select_centered_neighbor_pair:BadInput', ...
            'sys must be a struct.');
    end

    if ~isfield(sys, 'super_lattice') || isempty(sys.super_lattice)
        error('builder:select_centered_neighbor_pair:MissingSuperLattice', ...
            'sys.super_lattice is required and missing/empty.');
    end

    if ~isfield(sys, 'molecule_table') || isempty(sys.molecule_table)
        error('builder:select_centered_neighbor_pair:MissingMoleculeTable', ...
            'sys.molecule_table is required and missing/empty.');
    end
end