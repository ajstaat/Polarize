function ax = plot_local_environment(sys, focusMolIDs, varargin)
%VIZ.PLOT_LOCAL_ENVIRONMENT
% Plot focus molecules as WHOLE reconstructed molecule images from the
% original CONTCAR, while plotting neighbors faintly from sys for context.
%
% Required:
%   sys
%   focusMolIDs
%
% Name-value:
%   'Filename'             : path to CONTCAR/POSCAR (required for whole focus molecules)
%   'BondScale'            : default 1.20
%   'SortMolecules'        : default false
%   'NeighborShell'        : default 8.0
%   'ShowNeighbors'        : default true
%   'ShowBackground'       : default false
%   'ShowCOM'              : default true
%   'ShowMolIDs'           : default true
%   'ShowLabels'           : default false
%   'HighlightChargedSites': default true
%   'DrawCellBox'          : default false
%   'FocusMarkerSize'      : default 80
%   'NeighborMarkerSize'   : default 24
%   'BackgroundMarkerSize' : default 10

    p = inputParser;
    addRequired(p, 'sys', @isstruct);
    addRequired(p, 'focusMolIDs', @(x) isnumeric(x) && ~isempty(x));

    addParameter(p, 'Filename', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'BondScale', 1.20, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'SortMolecules', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'NeighborShell', 8.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'ShowNeighbors', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowBackground', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowCOM', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowMolIDs', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowLabels', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'HighlightChargedSites', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'DrawCellBox', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'FocusMarkerSize', 80, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'NeighborMarkerSize', 24, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'BackgroundMarkerSize', 10, @(x) isnumeric(x) && isscalar(x));
    parse(p, sys, focusMolIDs, varargin{:});
    opt = p.Results;

    focusMolIDs = unique(focusMolIDs(:));

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('sys.site_pos is missing.');
    end
    if ~isfield(sys, 'site_mol_id') || isempty(sys.site_mol_id)
        error('sys.site_mol_id is missing.');
    end

    T = builder.list_molecules(sys);

    if ~all(ismember(focusMolIDs, T.unique_mol_id))
        error('One or more focusMolIDs are not present in builder.list_molecules(sys).');
    end

    % Load whole base molecules from source file if provided
    useWholeFocus = strlength(string(opt.Filename)) > 0;
    if useWholeFocus
        mols = io.unwrap_all_contcar_molecules(opt.Filename, ...
            'BondScale', opt.BondScale, ...
            'SortMolecules', opt.SortMolecules);

        if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
            Hcell = sys.super_lattice;
        elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
            Hcell = sys.lattice;
        else
            error('Need sys.super_lattice or sys.lattice.');
        end

        % Recover the unit-cell lattice from the source file, not the supercell
        Ssrc = io.read_vasp_structure(opt.Filename);
        Hunit = Ssrc.lattice;
    end

    % Determine neighbors by COM
    neighborMolIDs = [];
    if opt.ShowNeighbors
        comAll = [T.cx, T.cy, T.cz];
        keepNeighbor = false(height(T), 1);

        for k = 1:numel(focusMolIDs)
            row = T(T.unique_mol_id == focusMolIDs(k), :);
            com0 = [row.cx, row.cy, row.cz];
            d = sqrt(sum((comAll - com0).^2, 2));
            keepNeighbor = keepNeighbor | (d <= opt.NeighborShell);
        end

        neighborMolIDs = T.unique_mol_id(keepNeighbor);
        neighborMolIDs = setdiff(neighborMolIDs, focusMolIDs);
    end

    backgroundMolIDs = setdiff(T.unique_mol_id, [focusMolIDs; neighborMolIDs]);

    figure('Color', 'w');
    ax = axes();
    hold(ax, 'on');

    if opt.DrawCellBox
        if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
            viz.draw_cell_box(ax, sys.super_lattice, [0 0 0], 'rows');
        elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
            viz.draw_cell_box(ax, sys.lattice, [0 0 0], 'rows');
        end
    end

    % Background omitted by default
    if opt.ShowBackground && ~isempty(backgroundMolIDs)
        idxBg = ismember(sys.site_mol_id, backgroundMolIDs);
        xyzBg = sys.site_pos(idxBg, :);
        scatter3(ax, xyzBg(:,1), xyzBg(:,2), xyzBg(:,3), ...
            opt.BackgroundMarkerSize, ...
            'MarkerFaceColor', [0.92 0.92 0.92], ...
            'MarkerEdgeColor', [0.92 0.92 0.92]);
    end

    % Neighbors faint gray
    if opt.ShowNeighbors && ~isempty(neighborMolIDs)
        for k = 1:numel(neighborMolIDs)
            idxN = (sys.site_mol_id == neighborMolIDs(k));
            xyzN = sys.site_pos(idxN, :);

            scatter3(ax, xyzN(:,1), xyzN(:,2), xyzN(:,3), ...
                opt.NeighborMarkerSize, ...
                'MarkerFaceColor', [0.85 0.85 0.85], ...
                'MarkerEdgeColor', [0.72 0.72 0.72]);
        end
    end

    % Focus molecules: reconstructed whole images if Filename is provided
    focusColors = lines(max(numel(focusMolIDs), 1));

    for k = 1:numel(focusMolIDs)
        uid = focusMolIDs(k);
        row = T(T.unique_mol_id == uid, :);

        if useWholeFocus
            baseID = row.base_mol_id;
            shift = [row.ix, row.iy, row.iz];

            Mol = mols{baseID};
            fracImg = Mol.fracUnwrapped + shift;
            xyzF = fracImg * Hunit;
        else
            idxF = (sys.site_mol_id == uid);
            xyzF = sys.site_pos(idxF, :);
        end

        scatter3(ax, xyzF(:,1), xyzF(:,2), xyzF(:,3), ...
            opt.FocusMarkerSize, ...
            'MarkerFaceColor', focusColors(k,:), ...
            'MarkerEdgeColor', 'k');

        com = mean(xyzF, 1);

        if opt.ShowCOM
            plot3(ax, com(1), com(2), com(3), 'kp', ...
                'MarkerSize', 12, 'MarkerFaceColor', 'y');
        end

        if opt.ShowMolIDs
            text(ax, com(1), com(2), com(3), sprintf('  %d', uid), ...
                'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k');
        end

        if opt.ShowLabels && useWholeFocus
            for i = 1:size(xyzF,1)
                text(ax, xyzF(i,1), xyzF(i,2), xyzF(i,3), ...
                    [' ' Mol.labels{i}], 'FontSize', 8);
            end
        end
    end

    % Highlight charged sites on focus molecules only.
    % These are plotted from sys.site_pos, so treat them as approximate locators.
    if opt.HighlightChargedSites && isfield(sys, 'site_charge') && ~isempty(sys.site_charge)
        idxFocus = ismember(sys.site_mol_id, focusMolIDs);
        idxCharged = abs(sys.site_charge) > 1e-14;
        idxQC = idxFocus & idxCharged;

        xyzQ = sys.site_pos(idxQC, :);
        q = sys.site_charge(idxQC);

        idxPos = q > 0;
        idxNeg = q < 0;

        if any(idxPos)
            scatter3(ax, xyzQ(idxPos,1), xyzQ(idxPos,2), xyzQ(idxPos,3), ...
                180, 'o', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
        end
        if any(idxNeg)
            scatter3(ax, xyzQ(idxNeg,1), xyzQ(idxNeg,2), xyzQ(idxNeg,3), ...
                180, 's', 'MarkerEdgeColor', 'b', 'LineWidth', 2);
        end
    end

    axis(ax, 'tight');
    grid(ax, 'on');
    xlabel(ax, 'x (A)');
    ylabel(ax, 'y (A)');
    zlabel(ax, 'z (A)');

    if numel(focusMolIDs) == 1
        title(ax, sprintf('Local environment around molecule %d', focusMolIDs));
    else
        title(ax, sprintf('Local environment around molecules %s', mat2str(focusMolIDs(:).')));
    end

    view(ax, 3);
    hold(ax, 'off');
end