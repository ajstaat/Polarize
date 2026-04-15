function ax = plot_molecule_pair(sys, uidA, uidB, varargin)
%VIZ.PLOT_MOLECULE_PAIR Plot only two unique molecule images.
% Molecules are unwrapped on the fly for visualization only.
%
% This does NOT modify sys or the calculation geometry.

    p = inputParser;
    addParameter(p, 'ShowLabels', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowCOM', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowMolIDs', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'HighlightChargedSites', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'MarkerSize', 70, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'ChargedMarkerSize', 180, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'BondScale', 1.20, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, varargin{:});
    opt = p.Results;

    if ~isfield(sys, 'site_pos') || isempty(sys.site_pos)
        error('sys.site_pos is missing.');
    end
    if ~isfield(sys, 'site_mol_id') || isempty(sys.site_mol_id)
        error('sys.site_mol_id is missing.');
    end
    if ~isfield(sys, 'site_type') || isempty(sys.site_type)
        error('sys.site_type is missing.');
    end

    H = get_plot_lattice(sys);

    idxA = find(sys.site_mol_id == uidA);
    idxB = find(sys.site_mol_id == uidB);

    if isempty(idxA), error('No sites found for uidA=%d', uidA); end
    if isempty(idxB), error('No sites found for uidB=%d', uidB); end

    [xyzA, orderA] = unwrap_molecule_for_plot(sys.site_pos(idxA,:), sys.site_type(idxA), H, opt.BondScale);
    [xyzB, orderB] = unwrap_molecule_for_plot(sys.site_pos(idxB,:), sys.site_type(idxB), H, opt.BondScale);

    idxA = idxA(orderA);
    idxB = idxB(orderB);

    figure('Color','w');
    ax = axes();
    hold(ax, 'on');

    scatter3(ax, xyzA(:,1), xyzA(:,2), xyzA(:,3), opt.MarkerSize, ...
        'MarkerFaceColor', [0.85 0.65 0.10], ...
        'MarkerEdgeColor', 'k');

    scatter3(ax, xyzB(:,1), xyzB(:,2), xyzB(:,3), opt.MarkerSize, ...
        'MarkerFaceColor', [0.20 0.45 0.85], ...
        'MarkerEdgeColor', 'k');

    comA = mean(xyzA, 1);
    comB = mean(xyzB, 1);

    if opt.ShowCOM
        plot3(ax, comA(1), comA(2), comA(3), 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
        plot3(ax, comB(1), comB(2), comB(3), 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'c');
        plot3(ax, [comA(1) comB(1)], [comA(2) comB(2)], [comA(3) comB(3)], ...
            'k--', 'LineWidth', 1.5);
    end

    if opt.ShowMolIDs
        text(ax, comA(1), comA(2), comA(3), sprintf('  %d', uidA), ...
            'FontSize', 11, 'FontWeight', 'bold');
        text(ax, comB(1), comB(2), comB(3), sprintf('  %d', uidB), ...
            'FontSize', 11, 'FontWeight', 'bold');
    end

    if opt.HighlightChargedSites && isfield(sys, 'site_charge') && ~isempty(sys.site_charge)
        qA = sys.site_charge(idxA);
        qB = sys.site_charge(idxB);

        idxPosA = qA > 1e-14;
        idxNegA = qA < -1e-14;
        idxPosB = qB > 1e-14;
        idxNegB = qB < -1e-14;

        if any(idxPosA)
            scatter3(ax, xyzA(idxPosA,1), xyzA(idxPosA,2), xyzA(idxPosA,3), ...
                opt.ChargedMarkerSize, 'o', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
        end
        if any(idxNegA)
            scatter3(ax, xyzA(idxNegA,1), xyzA(idxNegA,2), xyzA(idxNegA,3), ...
                opt.ChargedMarkerSize, 's', 'MarkerEdgeColor', 'b', 'LineWidth', 2);
        end
        if any(idxPosB)
            scatter3(ax, xyzB(idxPosB,1), xyzB(idxPosB,2), xyzB(idxPosB,3), ...
                opt.ChargedMarkerSize, 'o', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
        end
        if any(idxNegB)
            scatter3(ax, xyzB(idxNegB,1), xyzB(idxNegB,2), xyzB(idxNegB,3), ...
                opt.ChargedMarkerSize, 's', 'MarkerEdgeColor', 'b', 'LineWidth', 2);
        end
    end

    if opt.ShowLabels && isfield(sys, 'site_label') && ~isempty(sys.site_label)
        labelsA = sys.site_label(idxA);
        labelsB = sys.site_label(idxB);

        for i = 1:size(xyzA,1)
            text(ax, xyzA(i,1), xyzA(i,2), xyzA(i,3), [' ' char(labelsA{i})], 'FontSize', 8);
        end
        for i = 1:size(xyzB,1)
            text(ax, xyzB(i,1), xyzB(i,2), xyzB(i,3), [' ' char(labelsB{i})], 'FontSize', 8);
        end
    end

    axis(ax, 'equal');
    grid(ax, 'on');
    xlabel(ax, 'x (A)');
    ylabel(ax, 'y (A)');
    zlabel(ax, 'z (A)');
    title(ax, sprintf('Molecule pair: %d vs %d', uidA, uidB));
    view(ax, 3);
    hold(ax, 'off');
end

% -------------------------------------------------------------------------
function H = get_plot_lattice(sys)
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        H = sys.super_lattice;
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        H = sys.lattice;
    else
        error('Need sys.super_lattice or sys.lattice to unwrap molecule for plotting.');
    end

    if ~isequal(size(H), [3 3])
        error('Plot lattice must be 3x3.');
    end
end

% -------------------------------------------------------------------------
function [xyzUnwrapped, order] = unwrap_molecule_for_plot(xyzWrapped, types, H, bondScale)
% Unwrap one molecule for display using a PBC-aware bond graph and BFS.

    n = size(xyzWrapped, 1);
    order = (1:n).';

    fracWrapped = xyzWrapped / H;

    A = false(n, n);
    for i = 1:n-1
        ri = covalent_radius(types{i});
        for j = i+1:n
            rj = covalent_radius(types{j});

            df = fracWrapped(i,:) - fracWrapped(j,:);
            df = df - round(df);
            dcart = norm(df * H);

            cutoff = bondScale * (ri + rj);

            if strcmpi(types{i}, 'H') && strcmpi(types{j}, 'H')
                continue;
            end

            if dcart < cutoff
                A(i,j) = true;
                A(j,i) = true;
            end
        end
    end

    % If graph is disconnected, use connected component ordering and unwrap
    % each component separately, then shift all components near the first.
    G = graph(A);
    comp = conncomp(G);
    comps = unique(comp, 'stable');

    fracOut = nan(n,3);

    refCOM = [];
    for c = 1:numel(comps)
        idx = find(comp == comps(c));
        fracPart = unwrap_component(fracWrapped(idx,:), A(idx,idx));

        if isempty(refCOM)
            fracOut(idx,:) = fracPart;
            refCOM = mean(fracPart, 1);
        else
            partCOM = mean(fracPart, 1);
            d = partCOM - refCOM;
            d = d - round(d);
            shift = refCOM + d - partCOM;
            fracOut(idx,:) = fracPart + shift;
        end
    end

    xyzUnwrapped = fracOut * H;
end

% -------------------------------------------------------------------------
function fracUnwrapped = unwrap_component(fracWrapped, Aloc)
    n = size(fracWrapped, 1);

    fracUnwrapped = nan(n,3);
    visited = false(n,1);

    fracUnwrapped(1,:) = fracWrapped(1,:);
    visited(1) = true;
    queue = 1;

    while ~isempty(queue)
        current = queue(1);
        queue(1) = [];

        nbrs = find(Aloc(current,:));
        for nb = nbrs
            if ~visited(nb)
                deltaWrapped = fracWrapped(nb,:) - fracWrapped(current,:);
                deltaMin = deltaWrapped - round(deltaWrapped);
                fracUnwrapped(nb,:) = fracUnwrapped(current,:) + deltaMin;

                visited(nb) = true;
                queue(end+1) = nb; %#ok<AGROW>
            end
        end
    end

    for i = 1:n
        if ~visited(i)
            fracUnwrapped(i,:) = fracWrapped(i,:);
        end
    end
end

% -------------------------------------------------------------------------
function r = covalent_radius(sym)
    switch upper(sym)
        case 'H'
            r = 0.31;
        case 'C'
            r = 0.76;
        case 'N'
            r = 0.71;
        case 'O'
            r = 0.66;
        case 'F'
            r = 0.57;
        case 'S'
            r = 1.05;
        case 'CL'
            r = 1.02;
        case 'BR'
            r = 1.20;
        case 'I'
            r = 1.39;
        case 'SI'
            r = 1.11;
        otherwise
            error('No covalent radius defined for element: %s', sym);
    end
end