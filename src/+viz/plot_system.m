function ax = plot_system(sys, varargin)
%VIZ.PLOT_SYSTEM Plot atoms in the working Polarize system.
%
% Examples:
%   viz.plot_system(sys)
%   viz.plot_system(sys, 'ColorBy', 'unique_mol_id')
%   viz.plot_system(sys, 'ColorBy', 'base_mol_id', 'ShowCOM', true)
%   viz.plot_system(sys, 'ColorBy', 'unique_mol_id', ...
%       'ShowCOM', true, 'ShowMolIDs', true, 'DrawCellBox', true)

    p = inputParser;
    addParameter(p, 'ColorBy', 'unique_mol_id', @(x) ischar(x) || isstring(x));
    addParameter(p, 'MoleculeIDs', [], @(x) isnumeric(x));
    addParameter(p, 'ShowLabels', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'MarkerSize', 60, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'ShowCOM', false, @(x) islogical(x) && isscalar(x));

    % New options
    addParameter(p, 'ShowMolIDs', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'DrawCellBox', false, @(x) islogical(x) && isscalar(x));

    parse(p, varargin{:});
    opt = p.Results;

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('sys.site_pos is missing.');
    end

    xyz = sys.site_pos;

    switch char(opt.ColorBy)
        case 'unique_mol_id'
            if ~isfield(sys, 'site_mol_id') || isempty(sys.site_mol_id)
                error('sys.site_mol_id is missing.');
            end
            gid = sys.site_mol_id(:);
            titleStr = 'Colored by unique molecule ID';

        case 'base_mol_id'
            if ~isfield(sys, 'base_mol_id') || isempty(sys.base_mol_id)
                error('sys.base_mol_id is missing.');
            end
            gid = sys.base_mol_id(:);
            titleStr = 'Colored by base molecule ID';

        case 'site_type'
            if ~isfield(sys, 'site_type') || isempty(sys.site_type)
                error('sys.site_type is missing.');
            end
            [~,~,gid] = unique(string(sys.site_type));
            titleStr = 'Colored by site type';

        otherwise
            error('Unknown ColorBy option: %s', string(opt.ColorBy));
    end

    keep = true(size(gid));
    if ~isempty(opt.MoleculeIDs)
        keep = ismember(gid, opt.MoleculeIDs);
    end

    xyz = xyz(keep,:);
    gid = gid(keep);

    figure('Color','w');
    ax = axes();
    hold(ax, 'on');

    % Optional cell/supercell box
    if opt.DrawCellBox
        if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
            viz.draw_cell_box(ax, sys.super_lattice, [0 0 0]);
        elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
            viz.draw_cell_box(ax, sys.lattice, [0 0 0]);
        end
    end

    ug = unique(gid, 'stable');
    cmap = lines(max(numel(ug), 1));

    for k = 1:numel(ug)
        idx = (gid == ug(k));
        scatter3(ax, xyz(idx,1), xyz(idx,2), xyz(idx,3), opt.MarkerSize, ...
            'MarkerFaceColor', cmap(k,:), ...
            'MarkerEdgeColor', 'k', ...
            'DisplayName', sprintf('%s %d', char(opt.ColorBy), ug(k)));
    end

    if opt.ShowLabels
        if ~isfield(sys, 'site_label') || isempty(sys.site_label)
            warning('sys.site_label is missing; cannot show site labels.');
        else
            labels = sys.site_label(keep);
            for i = 1:size(xyz,1)
                text(ax, xyz(i,1), xyz(i,2), xyz(i,3), [' ' char(labels{i})], ...
                    'FontSize', 8);
            end
        end
    end

    if opt.ShowCOM || opt.ShowMolIDs
        T = builder.list_molecules(sys);

        % If plotting a subset, only show COMs/IDs for those molecules
        if ~isempty(opt.MoleculeIDs)
            switch char(opt.ColorBy)
                case {'unique_mol_id', 'base_mol_id'}
                    T = T(ismember(T.unique_mol_id, opt.MoleculeIDs) | ...
                          ismember(T.base_mol_id, opt.MoleculeIDs), :);
            end
        end

        if opt.ShowCOM
            for m = 1:height(T)
                plot3(ax, T.cx(m), T.cy(m), T.cz(m), 'kp', ...
                    'MarkerSize', 12, 'MarkerFaceColor', 'y');
            end
        end

        if opt.ShowMolIDs
            for m = 1:height(T)
                text(ax, T.cx(m), T.cy(m), T.cz(m), ...
                    sprintf(' %d', T.unique_mol_id(m)), ...
                    'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');
            end
        end
    end

    axis(ax, 'tight');
    daspect(ax, [1 1 1]);   % optional, try with/without
    grid(ax, 'on');
    xlabel(ax, 'x (A)');
    ylabel(ax, 'y (A)');
    zlabel(ax, 'z (A)');
    title(ax, titleStr);
    view(ax, 3);
    %legend(ax, 'Location', 'bestoutside');
    hold(ax, 'off');

    xlimVals = xlim(ax);
    ylimVals = ylim(ax);
    zlimVals = zlim(ax);

    fprintf('Plot bounds:\n');
    fprintf('  x: [%.4f, %.4f]\n', xlimVals(1), xlimVals(2));
    fprintf('  y: [%.4f, %.4f]\n', ylimVals(1), ylimVals(2));
    fprintf('  z: [%.4f, %.4f]\n', zlimVals(1), zlimVals(2));
end