function desc = complete_molecule_descriptors_relative_to_reference(sys, refMolID, varargin)
%COMPLETE_MOLECULE_DESCRIPTORS_RELATIVE_TO_REFERENCE
% Build a descriptor table for complete molecules relative to a reference.
%
% desc = builder.complete_molecule_descriptors_relative_to_reference(sys, refMolID)
% desc = builder.complete_molecule_descriptors_relative_to_reference(sys, refMolID, ...)
%
% Inputs
%   sys       working system struct from builder.make_crystal_system
%   refMolID  scalar supercell molecule ID
%
% Optional name-value inputs
%   'StackAxis'        : 'a' | 'b' | 'c' | numeric 1x3 vector, default 'b'
%   'IncludeReference' : logical, default false
%   'IncludeNormals'   : logical, default true
%   'Verbose'          : logical, default false
%
% Output
%   desc struct with fields:
%       .reference_mol_id
%       .reference_site_indices
%       .reference_com
%       .reference_normal
%       .stack_axis_hat
%       .table
%
% table columns
%   molecule_id
%   is_complete
%   n_sites
%   com_x, com_y, com_z
%   dr_x, dr_y, dr_z
%   d_par
%   d_perp_x, d_perp_y, d_perp_z
%   d_perp
%   distance
%   normal_x, normal_y, normal_z
%   normal_angle_deg

    p = inputParser;
    addRequired(p, 'sys', @isstruct);
    addRequired(p, 'refMolID', @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'StackAxis', 'b', @(x) (ischar(x) || isstring(x) || (isnumeric(x) && numel(x) == 3)));
    addParameter(p, 'IncludeReference', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'IncludeNormals', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
    parse(p, sys, refMolID, varargin{:});
    opt = p.Results;

    validate_sys(sys);
    validate_reference(sys, refMolID);

    T = sys.molecule_table;
    refRow = find(T.molecule_id == refMolID, 1, 'first');
    refIdx = builder.site_indices_for_molecule(sys, refMolID);
    refCOM = T.com(refRow, :);

    stackAxisHat = resolve_stack_axis(sys, opt.StackAxis);

    if opt.IncludeNormals
        refNormal = estimate_molecule_normal(sys.site_pos(refIdx, :));
    else
        refNormal = [NaN NaN NaN];
    end

    if opt.IncludeReference
        rows = (1:numel(T.molecule_id)).';
    else
        rows = find(T.molecule_id ~= refMolID);
    end

    rows = rows(T.is_complete_in_display(rows));
    nMol = numel(rows);

    molecule_id = zeros(nMol, 1);
    is_complete = false(nMol, 1);
    n_sites = zeros(nMol, 1);

    com = zeros(nMol, 3);
    dr = zeros(nMol, 3);
    d_par = zeros(nMol, 1);
    d_perp_vec = zeros(nMol, 3);
    d_perp = zeros(nMol, 1);
    distance = zeros(nMol, 1);

    normal = nan(nMol, 3);
    normal_angle_deg = nan(nMol, 1);

    for k = 1:nMol
        r = rows(k);
        thisMolID = T.molecule_id(r);
        idx = builder.site_indices_for_molecule(sys, thisMolID);

        thisCOM = T.com(r, :);
        thisDR = thisCOM - refCOM;
        thisDpar = dot(thisDR, stackAxisHat);
        thisDperpVec = thisDR - thisDpar * stackAxisHat;

        molecule_id(k) = thisMolID;
        is_complete(k) = T.is_complete_in_display(r);
        n_sites(k) = T.n_sites(r);

        com(k, :) = thisCOM;
        dr(k, :) = thisDR;
        d_par(k) = thisDpar;
        d_perp_vec(k, :) = thisDperpVec;
        d_perp(k) = norm(thisDperpVec);
        distance(k) = norm(thisDR);

        if opt.IncludeNormals
            thisNormal = estimate_molecule_normal(sys.site_pos(idx, :));
            normal(k, :) = thisNormal;
            normal_angle_deg(k) = unoriented_angle_deg(refNormal, thisNormal);
        end
    end

    [~, order] = sort(distance, 'ascend');

    desc = struct();
    desc.reference_mol_id = refMolID;
    desc.reference_site_indices = refIdx;
    desc.reference_com = refCOM;
    desc.reference_normal = refNormal;
    desc.stack_axis_hat = stackAxisHat;

    desc.table = table( ...
        molecule_id(order), ...
        is_complete(order), ...
        n_sites(order), ...
        com(order,1), com(order,2), com(order,3), ...
        dr(order,1), dr(order,2), dr(order,3), ...
        d_par(order), ...
        d_perp_vec(order,1), d_perp_vec(order,2), d_perp_vec(order,3), ...
        d_perp(order), ...
        distance(order), ...
        normal(order,1), normal(order,2), normal(order,3), ...
        normal_angle_deg(order), ...
        'VariableNames', { ...
            'molecule_id', ...
            'is_complete', ...
            'n_sites', ...
            'com_x', 'com_y', 'com_z', ...
            'dr_x', 'dr_y', 'dr_z', ...
            'd_par', ...
            'd_perp_x', 'd_perp_y', 'd_perp_z', ...
            'd_perp', ...
            'distance', ...
            'normal_x', 'normal_y', 'normal_z', ...
            'normal_angle_deg'});

    if opt.Verbose
        fprintf('Descriptor summary:\n');
        fprintf('  reference molecule ID   = %d\n', refMolID);
        fprintf('  complete molecules used = %d\n', nMol);
        fprintf('  stack axis hat          = [%9.4f %9.4f %9.4f]\n', stackAxisHat);
        if opt.IncludeNormals
            fprintf('  reference normal        = [%9.4f %9.4f %9.4f]\n', refNormal);
        end
    end
end


function validate_sys(sys)
    if ~isstruct(sys)
        error('builder:complete_molecule_descriptors_relative_to_reference:BadInput', ...
            'sys must be a struct.');
    end

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('builder:complete_molecule_descriptors_relative_to_reference:MissingSitePos', ...
            'sys.site_pos is required and missing/empty.');
    end

    if ~isfield(sys, 'molecule_table') || isempty(sys.molecule_table)
        error('builder:complete_molecule_descriptors_relative_to_reference:MissingMoleculeTable', ...
            'sys.molecule_table is required and missing/empty.');
    end

    T = sys.molecule_table;
    required = {'molecule_id', 'com', 'n_sites', 'is_complete_in_display'};
    for k = 1:numel(required)
        name = required{k};
        if ~isfield(T, name) || isempty(T.(name))
            error('builder:complete_molecule_descriptors_relative_to_reference:BadMoleculeTable', ...
                'sys.molecule_table.%s is required and missing/empty.', name);
        end
    end

    if size(T.com, 2) ~= 3
        error('builder:complete_molecule_descriptors_relative_to_reference:BadCOM', ...
            'sys.molecule_table.com must be N x 3.');
    end
end


function validate_reference(sys, refMolID)
    allMolIDs = sys.molecule_table.molecule_id(:);
    if ~ismember(refMolID, allMolIDs)
        error('builder:complete_molecule_descriptors_relative_to_reference:BadReferenceMolID', ...
            'Reference molecule ID %d is not present in sys.molecule_table.', refMolID);
    end
end


function stackAxisHat = resolve_stack_axis(sys, axisSpec)
    if ischar(axisSpec) || (isstring(axisSpec) && isscalar(axisSpec))
        axisSpec = lower(char(string(axisSpec)));

        if ~isfield(sys, 'super_lattice') || isempty(sys.super_lattice)
            error('builder:complete_molecule_descriptors_relative_to_reference:MissingSuperLattice', ...
                'sys.super_lattice is required for symbolic stack-axis selection.');
        end

        switch axisSpec
            case 'a'
                v = sys.super_lattice(1, :);
            case 'b'
                v = sys.super_lattice(2, :);
            case 'c'
                v = sys.super_lattice(3, :);
            otherwise
                error('builder:complete_molecule_descriptors_relative_to_reference:BadStackAxis', ...
                    'StackAxis must be ''a'', ''b'', ''c'', or a numeric 1x3 vector.');
        end
    elseif isnumeric(axisSpec) && numel(axisSpec) == 3
        v = reshape(axisSpec, 1, 3);
    else
        error('builder:complete_molecule_descriptors_relative_to_reference:BadStackAxis', ...
            'StackAxis must be ''a'', ''b'', ''c'', or a numeric 1x3 vector.');
    end

    nv = norm(v);
    if nv == 0
        error('builder:complete_molecule_descriptors_relative_to_reference:ZeroStackAxis', ...
            'Resolved stack axis has zero norm.');
    end

    stackAxisHat = v / nv;
end


function nHat = estimate_molecule_normal(X)
% Best-fit plane normal from site positions X (N x 3).
% Returns NaNs if the molecule does not provide enough geometric rank.

    if size(X, 1) < 3
        nHat = [NaN NaN NaN];
        return;
    end

    Xc = X - mean(X, 1);
    [~, S, V] = svd(Xc, 'econ');

    sing = diag(S);
    if numel(sing) < 3
        nHat = [NaN NaN NaN];
        return;
    end

    nHat = V(:, end).';
    nn = norm(nHat);
    if nn == 0
        nHat = [NaN NaN NaN];
    else
        nHat = nHat / nn;
    end
end


function ang = unoriented_angle_deg(n1, n2)
    if any(isnan(n1)) || any(isnan(n2))
        ang = NaN;
        return;
    end

    c = abs(dot(n1, n2));
    c = min(max(c, -1), 1);
    ang = acosd(c);
end