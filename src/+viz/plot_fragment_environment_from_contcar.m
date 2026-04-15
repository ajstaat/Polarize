function ax = plot_fragment_environment_from_contcar(filename, varargin)
%VIZ.PLOT_FRAGMENT_ENVIRONMENT_FROM_CONTCAR
% Plot a local molecular environment from whole reconstructed molecule
% images based on the original CONTCAR/POSCAR.
%
% This avoids mixed wrapped/unwrapped geometries.
%
% Required name-value:
%   'RefMolID'   : base molecule ID in the unit cell
%
% Main options:
%   'SecondMolID'    : optional second base molecule ID to also force inside box
%   'BondScale'      : default 1.20
%   'SortMolecules'  : default false
%   'Buffer'         : default 0.10 fractional units
%   'COMCutoff'      : default 15.0 A
%   'SearchImages'   : default [2 2 2]
%   'AlignRefCOMToOrigin' : default true
%   'ShowLabels'     : default false
%   'ShowCOM'        : default true
%   'ShowMolIDs'     : default true
%
% Color scheme:
%   - focus molecule: orange
%   - optional second focus molecule: blue
%   - neighbor atoms inside box: medium gray
%   - same neighbor atoms outside box: very light gray

    p = inputParser;
    addRequired(p, 'filename', @(x) ischar(x) || isstring(x));

    addParameter(p, 'RefMolID', [], @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'SecondMolID', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
    addParameter(p, 'BondScale', 1.20, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'SortMolecules', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'Buffer', 0.10, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'COMCutoff', 15.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'SearchImages', [2 2 2], @(x) isnumeric(x) && numel(x) == 3);
    addParameter(p, 'AlignRefCOMToOrigin', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowLabels', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowCOM', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowMolIDs', true, @(x) islogical(x) && isscalar(x));
    parse(p, filename, varargin{:});
    opt = p.Results;

    if isempty(opt.RefMolID)
        error('You must provide ''RefMolID''.');
    end

    S = io.read_vasp_structure(filename);
    H = S.lattice;

    mols = io.unwrap_all_contcar_molecules(filename, ...
        'BondScale', opt.BondScale, ...
        'SortMolecules', opt.SortMolecules);

    nBase = numel(mols);

    if opt.RefMolID > nBase
        error('RefMolID exceeds number of molecules found (%d).', nBase);
    end
    if ~isempty(opt.SecondMolID) && opt.SecondMolID > nBase
        error('SecondMolID exceeds number of molecules found (%d).', nBase);
    end

    % ------------------------------------------------------------
    % Build enclosing box from the reference molecule (and optional second)
    % ------------------------------------------------------------
    fRef = mols{opt.RefMolID}.fracUnwrapped;

    if isempty(opt.SecondMolID)
        fUnion = fRef;
    else
        f2 = mols{opt.SecondMolID}.fracUnwrapped;
        fUnion = [fRef; f2];
    end

    fmin0 = min(fUnion, [], 1);
    shift = opt.Buffer - fmin0;

    fRefShift = fRef + shift;

    if isempty(opt.SecondMolID)
        fUnionShift = fRefShift;
    else
        f2Shift = mols{opt.SecondMolID}.fracUnwrapped + shift;
        fUnionShift = [fRefShift; f2Shift];
    end

    boxDims = ceil(max(fUnionShift, [], 1));

    % Cartesian box corners for plotting
    boxOriginFrac = [0 0 0];
    boxOriginCart = boxOriginFrac * H;
    boxLattice = diag(boxDims) * H;  % rows = box vectors in Cartesian

    % ------------------------------------------------------------
    % Generate nearby molecule images from integer lattice shifts
    % ------------------------------------------------------------
    sx = -opt.SearchImages(1):opt.SearchImages(1);
    sy = -opt.SearchImages(2):opt.SearchImages(2);
    sz = -opt.SearchImages(3):opt.SearchImages(3);

    refCOM = mean(fRefShift * H, 1);

    images = struct( ...
        'base_mol_id', {}, ...
        'shift', {}, ...
        'frac', {}, ...
        'cart', {}, ...
        'com', {}, ...
        'insideMask', {} );

    count = 0;
    for b = 1:nBase
        f0 = mols{b}.fracUnwrapped + shift;

        for ix = sx
            for iy = sy
                for iz = sz
                    sh = [ix iy iz];
                    fImg = f0 + sh;
                    cartImg = fImg * H;
                    comImg = mean(cartImg, 1);

                    if norm(comImg - refCOM) <= opt.COMCutoff
                        count = count + 1;
                        images(count).base_mol_id = b; %#ok<AGROW>
                        images(count).shift = sh;
                        images(count).frac = fImg;
                        images(count).cart = cartImg;
                        images(count).com = comImg;

                        inside = ...
                            fImg(:,1) >= 0 & fImg(:,1) <= boxDims(1) & ...
                            fImg(:,2) >= 0 & fImg(:,2) <= boxDims(2) & ...
                            fImg(:,3) >= 0 & fImg(:,3) <= boxDims(3);
                        images(count).insideMask = inside;
                    end
                end
            end
        end
    end

    % ------------------------------------------------------------
    % Optional global translation so ref COM is at origin for easier viewing
    % ------------------------------------------------------------
    globalShift = [0 0 0];
    if opt.AlignRefCOMToOrigin
        globalShift = -refCOM;
    end

    % ------------------------------------------------------------
    % Plot
    % ------------------------------------------------------------
    figure('Color', 'w');
    ax = axes();
    hold(ax, 'on');

    % Draw the enclosing box
    viz.draw_cell_box(ax, boxLattice, boxOriginCart + globalShift, 'rows');

    focusColor1 = [0.85 0.65 0.10];
    focusColor2 = [0.20 0.45 0.85];
    grayInside = [0.70 0.70 0.70];
    grayOutside = [0.90 0.90 0.90];

    for k = 1:numel(images)
        b = images(k).base_mol_id;
        xyz = images(k).cart + globalShift;
        com = images(k).com + globalShift;
        inside = images(k).insideMask;

        isRef = (b == opt.RefMolID);
        isSecond = ~isempty(opt.SecondMolID) && (b == opt.SecondMolID);

        if isRef
            scatter3(ax, xyz(:,1), xyz(:,2), xyz(:,3), 80, ...
                'MarkerFaceColor', focusColor1, ...
                'MarkerEdgeColor', 'k');

        elseif isSecond
            scatter3(ax, xyz(:,1), xyz(:,2), xyz(:,3), 80, ...
                'MarkerFaceColor', focusColor2, ...
                'MarkerEdgeColor', 'k');

        else
            if any(inside)
                scatter3(ax, xyz(inside,1), xyz(inside,2), xyz(inside,3), 22, ...
                    'MarkerFaceColor', grayInside, ...
                    'MarkerEdgeColor', grayInside);
            end
            if any(~inside)
                scatter3(ax, xyz(~inside,1), xyz(~inside,2), xyz(~inside,3), 16, ...
                    'MarkerFaceColor', grayOutside, ...
                    'MarkerEdgeColor', grayOutside);
            end
        end

        if opt.ShowCOM
            plot3(ax, com(1), com(2), com(3), 'kp', ...
                'MarkerSize', 8, 'MarkerFaceColor', 'y');
        end

        if opt.ShowMolIDs
            label = sprintf('b%d [%d %d %d]', ...
                images(k).base_mol_id, images(k).shift(1), images(k).shift(2), images(k).shift(3));
            text(ax, com(1), com(2), com(3), [' ' label], ...
                'FontSize', 8, 'Color', [0.2 0.2 0.2]);
        end
    end

    if opt.ShowLabels
        for k = 1:numel(images)
            b = images(k).base_mol_id;
            xyz = images(k).cart + globalShift;
            labels = mols{b}.labels;

            if b == opt.RefMolID || (~isempty(opt.SecondMolID) && b == opt.SecondMolID)
                for i = 1:size(xyz,1)
                    text(ax, xyz(i,1), xyz(i,2), xyz(i,3), [' ' labels{i}], ...
                        'FontSize', 8);
                end
            end
        end
    end

    axis(ax, 'tight');
    grid(ax, 'on');
    xlabel(ax, 'x (A)');
    ylabel(ax, 'y (A)');
    zlabel(ax, 'z (A)');
    title(ax, sprintf('Fragment environment around base molecule %d', opt.RefMolID));
    view(ax, 3);
    hold(ax, 'off');
end