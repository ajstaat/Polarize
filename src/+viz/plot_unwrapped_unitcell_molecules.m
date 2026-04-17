function out = plot_unwrapped_unitcell_molecules(filename, varargin)
%PLOT_UNWRAPPED_UNITCELL_MOLECULES Plot all unwrapped molecules from one
% POSCAR/CONTCAR-style unit cell.
%
% out = viz.plot_unwrapped_unitcell_molecules(filename)
% out = viz.plot_unwrapped_unitcell_molecules(filename, ...)
%
% Optional name-value inputs
%   'BondScale'       : default 1.20
%   'SortMolecules'   : default false
%   'DrawBox'         : default true
%   'ShowCOM'         : default true
%   'ShowLabels'      : default true
%   'MarkerSize'      : default 72
%   'Axes'            : axes handle, default []
%   'Title'           : char/string, default auto
%
% Output
%   out struct with fields:
%       .ax
%       .molecules
%       .com
%       .lattice

    p = inputParser;
    addRequired(p, 'filename', @(x) ischar(x) || isstring(x));

    addParameter(p, 'BondScale', 1.20, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'SortMolecules', false, @(x) islogical(x) && isscalar(x));

    addParameter(p, 'DrawBox', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowCOM', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowLabels', true, @(x) islogical(x) && isscalar(x));

    addParameter(p, 'MarkerSize', 72, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Axes', [], @(x) isempty(x) || isgraphics(x, 'axes'));
    addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));

    parse(p, filename, varargin{:});
    opt = p.Results;

    S = io.read_vasp_structure(filename);
    mols = io.unwrap_all_contcar_molecules(filename, ...
        'BondScale', opt.BondScale, ...
        'SortMolecules', opt.SortMolecules);

    if isempty(opt.Axes)
        figure('Color', 'w');
        ax = axes();
    else
        ax = opt.Axes;
        cla(ax);
    end

    hold(ax, 'on');

    % Softer palette, closer to recent plot styling
    baseColors = [
        0.15 0.45 0.75
        0.90 0.55 0.15
        0.35 0.65 0.35
        0.55 0.35 0.75
        0.80 0.30 0.30
        0.20 0.65 0.65
        0.65 0.50 0.20
        0.45 0.45 0.45];

    nMol = numel(mols);
    com = zeros(nMol, 3);

    for k = 1:nMol
        xyz = mols{k}.cart;
        com(k, :) = mean(xyz, 1);

        c = baseColors(mod(k-1, size(baseColors,1)) + 1, :);

        scatter3(ax, ...
            xyz(:,1), xyz(:,2), xyz(:,3), ...
            opt.MarkerSize, ...
            c, ...
            'filled', ...
            'MarkerEdgeColor', [0.15 0.15 0.15]);

        if opt.ShowCOM
            scatter3(ax, ...
                com(k,1), com(k,2), com(k,3), ...
                95, ...
                c, ...
                'd', ...
                'filled', ...
                'MarkerEdgeColor', 'k');
        end

        if opt.ShowLabels
            offset = 0.012 * sum(abs(S.lattice), 1);

            text(ax, ...
                com(k,1) + offset(1), ...
                com(k,2) + offset(2), ...
                com(k,3) + offset(3), ...
                sprintf('mol %d', k), ...
                'FontWeight', 'bold', ...
                'Color', c);
        end
    end

    if opt.DrawBox
        viz.draw_cell_box(ax, S.lattice, [0 0 0], 'rows');
    end

    xlabel(ax, 'x');
    ylabel(ax, 'y');
    zlabel(ax, 'z');
    axis(ax, 'equal');
    grid(ax, 'on');
    view(ax, 3);

    if strlength(string(opt.Title)) > 0
        title(ax, char(string(opt.Title)));
    else
        title(ax, sprintf('All unwrapped molecules in unit cell (n = %d)', nMol));
    end

    hold(ax, 'off');

    out = struct();
    out.ax = ax;
    out.molecules = mols;
    out.com = com;
    out.lattice = S.lattice;
end