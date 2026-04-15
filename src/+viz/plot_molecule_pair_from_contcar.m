function ax = plot_molecule_pair_from_contcar(filename, varargin)
%VIZ.PLOT_MOLECULE_PAIR_FROM_CONTCAR
% General pair plotter for molecules from a CONTCAR/POSCAR.
%
% This function subsumes both:
%   - plot_contcar_pair
%   - plot_selected_pair_from_contcar
%
% It supports either:
%   (A) plotting two base molecules from the unit cell, or
%   (B) plotting a reference base molecule plus a selected translated image
%       from a descriptor/selector row.
%
% -------------------------------------------------------------------------
% Supported usage patterns
% -------------------------------------------------------------------------
%
% 1) Two base molecules in the unit cell:
%
%   viz.plot_molecule_pair_from_contcar(filename, ...
%       'RefMolID', 1, ...
%       'NeighborMolID', 2)
%
% This plots:
%   ref molecule      = base molecule 1 at [0 0 0]
%   neighbor molecule = base molecule 2 at [0 0 0]
%
% 2) Reference + selected neighbor image from selector output:
%
%   viz.plot_molecule_pair_from_contcar(filename, ...
%       'RefMolID', 1, ...
%       'NeighborRow', rowNeighbor)
%
% where rowNeighbor is a one-row table containing at least:
%   base_mol_id, ix, iy, iz
%
% This plots:
%   ref molecule      = base molecule 1 at RefShift
%   neighbor molecule = base_mol_id at [ix iy iz]
%
% -------------------------------------------------------------------------
% Required input
% -------------------------------------------------------------------------
%   filename
%
% Required name-value
%   'RefMolID' : base molecule ID of the reference molecule
%
% Optional name-value
%   'RefShift'        : default [0 0 0]
%   'NeighborMolID'   : base molecule ID for simple unit-cell pair mode
%   'NeighborShift'   : default [0 0 0]
%   'NeighborRow'     : one-row table with base_mol_id, ix, iy, iz
%   'BondScale'       : default 1.20
%   'SortMolecules'   : default false
%   'CenterMode'      : 'reference_com' | 'midpoint' | 'none'
%                       default = 'reference_com'
%   'ShowCOM'         : default true
%   'ShowLabels'      : default false
%   'ShowPairLine'    : default true
%   'ShowTitle'       : default true
%
% Notes
%   - If NeighborRow is provided, it takes precedence over NeighborMolID.
%   - This function uses whole unwrapped molecules from
%     io.unwrap_all_contcar_molecules, so the pair geometry is consistent.
%
% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
%   ax : axes handle

    p = inputParser;
    addRequired(p, 'filename', @(x) ischar(x) || isstring(x));

    addParameter(p, 'RefMolID', [], @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'RefShift', [0 0 0], @(x) isnumeric(x) && numel(x) == 3);

    addParameter(p, 'NeighborMolID', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
    addParameter(p, 'NeighborShift', [0 0 0], @(x) isnumeric(x) && numel(x) == 3);
    addParameter(p, 'NeighborRow', [], @(x) isempty(x) || (istable(x) && height(x) == 1));

    addParameter(p, 'BondScale', 1.20, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'SortMolecules', false, @(x) islogical(x) && isscalar(x));

    addParameter(p, 'CenterMode', 'reference_com', @(x) ischar(x) || isstring(x));
    addParameter(p, 'ShowCOM', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowLabels', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowPairLine', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowTitle', true, @(x) islogical(x) && isscalar(x));
    parse(p, filename, varargin{:});
    opt = p.Results;

    if isempty(opt.RefMolID)
        error('You must provide ''RefMolID''.');
    end

    % ---------------------------------------------------------------------
    % Load source structure and whole unwrapped base molecules
    % ---------------------------------------------------------------------
    S = io.read_vasp_structure(filename);
    H = S.lattice;

    mols = io.unwrap_all_contcar_molecules(filename, ...
        'BondScale', opt.BondScale, ...
        'SortMolecules', opt.SortMolecules);

    nBase = numel(mols);
    if opt.RefMolID > nBase
        error('RefMolID exceeds number of molecules found (%d).', nBase);
    end

    % ---------------------------------------------------------------------
    % Resolve neighbor specification
    % ---------------------------------------------------------------------
    refMolID = opt.RefMolID;
    refShift = reshape(opt.RefShift, 1, 3);

    if ~isempty(opt.NeighborRow)
        row = opt.NeighborRow;

        requiredCols = {'base_mol_id','ix','iy','iz'};
        missing = setdiff(requiredCols, row.Properties.VariableNames);
        if ~isempty(missing)
            error('NeighborRow is missing required columns: %s', strjoin(missing, ', '));
        end

        neighborMolID = row.base_mol_id(1);
        neighborShift = [row.ix(1), row.iy(1), row.iz(1)];
        neighborLabel = sprintf('b%d [%d %d %d]', ...
            neighborMolID, neighborShift(1), neighborShift(2), neighborShift(3));
    else
        if isempty(opt.NeighborMolID)
            error('Provide either ''NeighborMolID'' or ''NeighborRow''.');
        end

        neighborMolID = opt.NeighborMolID;
        neighborShift = reshape(opt.NeighborShift, 1, 3);
        neighborLabel = sprintf('b%d [%d %d %d]', ...
            neighborMolID, neighborShift(1), neighborShift(2), neighborShift(3));
    end

    if neighborMolID > nBase
        error('NeighborMolID exceeds number of molecules found (%d).', nBase);
    end

    % ---------------------------------------------------------------------
    % Reconstruct whole molecule images
    % ---------------------------------------------------------------------
    MolRef = mols{refMolID};
    fracRef = MolRef.fracUnwrapped + refShift;
    xyzRef = fracRef * H;
    comRef = mean(xyzRef, 1);

    MolNbr = mols{neighborMolID};
    fracNbr = MolNbr.fracUnwrapped + neighborShift;
    xyzNbr = fracNbr * H;
    comNbr = mean(xyzNbr, 1);

    % ---------------------------------------------------------------------
    % Recenter for display
    % ---------------------------------------------------------------------
    switch lower(string(opt.CenterMode))
        case "reference_com"
            shift = -comRef;

        case "midpoint"
            shift = -0.5 * (comRef + comNbr);

        case "none"
            shift = [0 0 0];

        otherwise
            error('CenterMode must be ''reference_com'', ''midpoint'', or ''none''.');
    end

    xyzRef = xyzRef + shift;
    xyzNbr = xyzNbr + shift;
    comRef = comRef + shift;
    comNbr = comNbr + shift;

    % ---------------------------------------------------------------------
    % Plot
    % ---------------------------------------------------------------------
    figure('Color', 'w');
    ax = axes();
    hold(ax, 'on');

    colorRef = [0.85 0.65 0.10];
    colorNbr = [0.20 0.45 0.85];

    scatter3(ax, xyzRef(:,1), xyzRef(:,2), xyzRef(:,3), 80, ...
        'MarkerFaceColor', colorRef, ...
        'MarkerEdgeColor', 'k');

    scatter3(ax, xyzNbr(:,1), xyzNbr(:,2), xyzNbr(:,3), 80, ...
        'MarkerFaceColor', colorNbr, ...
        'MarkerEdgeColor', 'k');

    if opt.ShowCOM
        plot3(ax, comRef(1), comRef(2), comRef(3), 'kp', ...
            'MarkerSize', 10, 'MarkerFaceColor', 'y');
        plot3(ax, comNbr(1), comNbr(2), comNbr(3), 'kp', ...
            'MarkerSize', 10, 'MarkerFaceColor', 'c');
    end

    if opt.ShowPairLine
        plot3(ax, [comRef(1) comNbr(1)], [comRef(2) comNbr(2)], [comRef(3) comNbr(3)], ...
            'k--', 'LineWidth', 1.3);
    end

    % Labels at COM
    text(ax, comRef(1), comRef(2), comRef(3), ...
        sprintf('  ref b%d [%d %d %d]', refMolID, refShift(1), refShift(2), refShift(3)), ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');

    text(ax, comNbr(1), comNbr(2), comNbr(3), ['  ' neighborLabel], ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');

    % Optional atom labels
    if opt.ShowLabels
        labelsRef = MolRef.labels;
        labelsNbr = MolNbr.labels;

        for i = 1:size(xyzRef,1)
            text(ax, xyzRef(i,1), xyzRef(i,2), xyzRef(i,3), [' ' labelsRef{i}], ...
                'FontSize', 8);
        end

        for i = 1:size(xyzNbr,1)
            text(ax, xyzNbr(i,1), xyzNbr(i,2), xyzNbr(i,3), [' ' labelsNbr{i}], ...
                'FontSize', 8);
        end
    end

    axis(ax, 'tight');
    grid(ax, 'on');
    xlabel(ax, 'x (A)');
    ylabel(ax, 'y (A)');
    zlabel(ax, 'z (A)');

    if opt.ShowTitle
        title(ax, 'Reference molecule + selected neighbor image');
    end

    view(ax, 3);
    hold(ax, 'off');
end