function ax = plot_all_unwrapped_contcar_molecules(filename, varargin)
%PLOT_ALL_UNWRAPPED_CONTCAR_MOLECULES Plot every unwrapped molecule from a CONTCAR.

    p = inputParser;
    addRequired(p, 'filename', @(x) ischar(x) || isstring(x));
    addParameter(p, 'BondScale', 1.20, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'SortMolecules', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowCOM', true, @(x) islogical(x) && isscalar(x));
    parse(p, filename, varargin{:});
    opt = p.Results;

    molecules = io.unwrap_all_contcar_molecules(filename, ...
        'BondScale', opt.BondScale, ...
        'SortMolecules', opt.SortMolecules);

    figure('Color','w');
    ax = axes();
    hold(ax, 'on');

    cmap = lines(max(numel(molecules), 1));

    for k = 1:numel(molecules)
        xyz = molecules{k}.cart;
        scatter3(ax, xyz(:,1), xyz(:,2), xyz(:,3), 55, ...
            'MarkerFaceColor', cmap(k,:), ...
            'MarkerEdgeColor', 'k');

        if opt.ShowCOM
            com = molecules{k}.com;
            plot3(ax, com(1), com(2), com(3), 'kp', ...
                'MarkerSize', 10, 'MarkerFaceColor', 'y');
            text(ax, com(1), com(2), com(3), sprintf('  %d', k), ...
                'FontSize', 10, 'FontWeight', 'bold');
        end
    end

    axis(ax, 'equal');
    grid(ax, 'on');
    xlabel(ax, 'x (A)');
    ylabel(ax, 'y (A)');
    zlabel(ax, 'z (A)');
    title(ax, 'All unwrapped CONTCAR molecules');
    view(ax, 3);
    hold(ax, 'off');
end