function out = plot_supercell_selection(sys, varargin)
%PLOT_SUPERCELL_SELECTION Plot a built supercell and highlight selected
% complete molecule(s) from the actual stored supercell coordinates.
%
% out = viz.plot_supercell_selection(sys)
% out = viz.plot_supercell_selection(sys, ...)
%
% Required input
%   sys   working system struct from builder.make_crystal_system
%
% Optional name-value inputs
%   'ReferenceMolID'        scalar supercell molecule ID, default []
%                           if empty, choose closest COMPLETE molecule to center
%   'NeighborMolID'         scalar supercell molecule ID, default []
%   'RequireCompleteRef'    logical, default true
%   'RequireCompleteNbr'    logical, default true
%   'IncompleteAction'      'error' | 'warning', default 'error'
%
%   'DrawBox'               logical, default true
%   'ShowCOM'               logical, default true
%   'ShowLabels'            logical, default true
%   'ShowPairLine'          logical, default true
%
%   'EnvMarkerSize'         default 52
%   'RefMarkerSize'         default 78
%   'NbrMarkerSize'         default 78
%
%   'EnvColor'              default [0.80 0.80 0.80]
%   'RefColor'              default [0.15 0.45 0.75]
%   'NbrColor'              default [0.90 0.55 0.15]
%
%   'EnvFaceAlpha'          default 0.16
%   'EnvEdgeAlpha'          default 0.08
%
%   'Axes'                  axes handle, default []
%   'Title'                 char/string, default auto
%
% Output
%   out struct with fields:
%       .ax
%       .reference_unique_mol_id
%       .neighbor_unique_mol_id
%       .reference_site_indices
%       .neighbor_site_indices
%       .reference_com
%       .neighbor_com
%       .supercell_center

    p = inputParser;
    addRequired(p, 'sys', @isstruct);

    addParameter(p, 'ReferenceMolID', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'NeighborMolID',  [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));

    addParameter(p, 'RequireCompleteRef', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'RequireCompleteNbr', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'IncompleteAction', 'error', @(x) ischar(x) || isstring(x));

    addParameter(p, 'DrawBox', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowCOM', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowLabels', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ShowPairLine', true, @(x) islogical(x) && isscalar(x));

    addParameter(p, 'EnvMarkerSize', 52, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'RefMarkerSize', 78, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'NbrMarkerSize', 78, @(x) isnumeric(x) && isscalar(x) && x > 0);

    addParameter(p, 'EnvColor', [0.80 0.80 0.80], @(x) isnumeric(x) && numel(x) == 3);
    addParameter(p, 'RefColor', [0.15 0.45 0.75], @(x) isnumeric(x) && numel(x) == 3);
    addParameter(p, 'NbrColor', [0.90 0.55 0.15], @(x) isnumeric(x) && numel(x) == 3);

    addParameter(p, 'EnvFaceAlpha', 0.16, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    addParameter(p, 'EnvEdgeAlpha', 0.08, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);

    addParameter(p, 'Axes', [], @(x) isempty(x) || isgraphics(x, 'axes'));
    addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));

    parse(p, sys, varargin{:});
    opt = p.Results;

    validate_sys(sys);

    incompleteAction = lower(char(string(opt.IncompleteAction)));
    if ~ismember(incompleteAction, {'error', 'warning'})
        error('viz:plot_supercell_selection:BadIncompleteAction', ...
            'IncompleteAction must be ''error'' or ''warning''.');
    end

    if isempty(opt.Axes)
        figure('Color', 'w');
        ax = axes();
    else
        ax = opt.Axes;
        cla(ax);
    end

    hold(ax, 'on');

    % ---------------------------------------------------------------------
    % Resolve reference molecule
    % ---------------------------------------------------------------------
    if isempty(opt.ReferenceMolID)
        refUID = choose_center_complete_molecule(sys);
    else
        refUID = double(opt.ReferenceMolID);
        validate_molecule_id(sys, refUID, 'ReferenceMolID');
    end

    refRow = molecule_table_row(sys.molecule_table, refUID);
    enforce_completeness_policy(refRow, 'reference', incompleteAction, opt.RequireCompleteRef);

    refIdx = builder.site_indices_for_molecule(sys, refUID);
    refCOM = mean(sys.site_pos(refIdx, :), 1);

    % ---------------------------------------------------------------------
    % Resolve neighbor molecule
    % ---------------------------------------------------------------------
    nbrUID = [];
    nbrIdx = [];
    nbrCOM = [];

    if ~isempty(opt.NeighborMolID)
        nbrUID = double(opt.NeighborMolID);
        validate_molecule_id(sys, nbrUID, 'NeighborMolID');

        if nbrUID == refUID
            error('viz:plot_supercell_selection:DuplicateSelection', ...
                'ReferenceMolID and NeighborMolID must be different.');
        end

        nbrRow = molecule_table_row(sys.molecule_table, nbrUID);
        enforce_completeness_policy(nbrRow, 'neighbor', incompleteAction, opt.RequireCompleteNbr);

        nbrIdx = builder.site_indices_for_molecule(sys, nbrUID);
        nbrCOM = mean(sys.site_pos(nbrIdx, :), 1);
    end

    % ---------------------------------------------------------------------
    % Plot environment cloud
    % ---------------------------------------------------------------------
    scatter3(ax, ...
        sys.site_pos(:,1), sys.site_pos(:,2), sys.site_pos(:,3), ...
        opt.EnvMarkerSize, opt.EnvColor, 'filled', ...
        'MarkerFaceAlpha', opt.EnvFaceAlpha, ...
        'MarkerEdgeAlpha', opt.EnvEdgeAlpha);

    % ---------------------------------------------------------------------
    % Plot selected molecule(s) from actual stored supercell coordinates
    % ---------------------------------------------------------------------
    scatter3(ax, ...
        sys.site_pos(refIdx,1), sys.site_pos(refIdx,2), sys.site_pos(refIdx,3), ...
        opt.RefMarkerSize, opt.RefColor, 'filled', ...
        'MarkerEdgeColor', [0.15 0.15 0.15]);

    if ~isempty(nbrUID)
        scatter3(ax, ...
            sys.site_pos(nbrIdx,1), sys.site_pos(nbrIdx,2), sys.site_pos(nbrIdx,3), ...
            opt.NbrMarkerSize, opt.NbrColor, 'filled', ...
            'MarkerEdgeColor', [0.15 0.15 0.15]);
    end

    % ---------------------------------------------------------------------
    % COM markers and optional connection line
    % ---------------------------------------------------------------------
    if opt.ShowCOM
        scatter3(ax, refCOM(1), refCOM(2), refCOM(3), ...
            90, opt.RefColor, 'd', 'filled', 'MarkerEdgeColor', 'k');

        if ~isempty(nbrUID)
            scatter3(ax, nbrCOM(1), nbrCOM(2), nbrCOM(3), ...
                90, opt.NbrColor, 's', 'filled', 'MarkerEdgeColor', 'k');
        end
    end

    if opt.ShowPairLine && ~isempty(nbrUID)
        plot3(ax, ...
            [refCOM(1), nbrCOM(1)], ...
            [refCOM(2), nbrCOM(2)], ...
            [refCOM(3), nbrCOM(3)], ...
            'k--', 'LineWidth', 1.1);
    end

    % ---------------------------------------------------------------------
    % Labels
    % ---------------------------------------------------------------------
    if opt.ShowLabels
        offset = 0.012 * sum(abs(sys.super_lattice), 1);

        text(ax, refCOM(1)+offset(1), refCOM(2)+offset(2), refCOM(3)+offset(3), ...
            sprintf('ref %d', refUID), ...
            'FontWeight', 'bold', 'Color', opt.RefColor);

        if ~isempty(nbrUID)
            text(ax, nbrCOM(1)+offset(1), nbrCOM(2)+offset(2), nbrCOM(3)+offset(3), ...
                sprintf('nbr %d', nbrUID), ...
                'FontWeight', 'bold', 'Color', opt.NbrColor);
        end
    end

    if opt.DrawBox
        viz.draw_cell_box(ax, sys.super_lattice, [0 0 0], 'rows');
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
        if isempty(nbrUID)
            title(ax, sprintf('Supercell selection: complete ref %d', refUID));
        else
            title(ax, sprintf('Supercell selection: ref %d, nbr %d', refUID, nbrUID));
        end
    end

    hold(ax, 'off');

    out = struct();
    out.ax = ax;
    out.reference_unique_mol_id = refUID;
    out.neighbor_unique_mol_id = nbrUID;
    out.reference_site_indices = refIdx;
    out.neighbor_site_indices = nbrIdx;
    out.reference_com = refCOM;
    out.neighbor_com = nbrCOM;
    out.supercell_center = 0.5 * sum(sys.super_lattice, 1);
end


function validate_sys(sys)
    required = {'site_pos', 'site_mol_id', 'super_lattice', 'molecule_table'};
    for k = 1:numel(required)
        name = required{k};
        if ~isfield(sys, name) || isempty(sys.(name))
            error('viz:plot_supercell_selection:MissingField', ...
                'sys.%s is required and missing/empty.', name);
        end
    end

    if size(sys.site_pos, 2) ~= 3
        error('viz:plot_supercell_selection:BadSitePos', ...
            'sys.site_pos must be N x 3.');
    end

    T = sys.molecule_table;
    reqTableFields = {'molecule_id', 'site_indices', 'n_sites', 'com', ...
                      'is_complete_in_display', 'n_display_fragments', ...
                      'largest_fragment_fraction'};
    for k = 1:numel(reqTableFields)
        name = reqTableFields{k};
        if ~isfield(T, name) || isempty(T.(name))
            error('viz:plot_supercell_selection:BadMoleculeTable', ...
                'sys.molecule_table.%s is required and missing/empty.', name);
        end
    end
end


function validate_molecule_id(sys, molID, argName)
    allMolIDs = unique(sys.site_mol_id(:)).';
    if ~ismember(molID, allMolIDs)
        error('viz:plot_supercell_selection:BadMoleculeID', ...
            '%s = %d is not present in sys.site_mol_id.', argName, molID);
    end
end


function row = molecule_table_row(T, molID)
    idx = find(T.molecule_id == molID, 1, 'first');
    if isempty(idx)
        error('viz:plot_supercell_selection:BadMoleculeID', ...
            'Molecule ID %d not found in sys.molecule_table.', molID);
    end

    row = struct();
    row.molecule_id = T.molecule_id(idx);
    row.site_indices = T.site_indices{idx};
    row.n_sites = T.n_sites(idx);
    row.com = T.com(idx, :);
    row.is_complete_in_display = T.is_complete_in_display(idx);
    row.n_display_fragments = T.n_display_fragments(idx);
    row.largest_fragment_fraction = T.largest_fragment_fraction(idx);
end


function enforce_completeness_policy(row, label, action, requireComplete)
    if ~requireComplete || row.is_complete_in_display
        return;
    end

    msg = sprintf([ ...
        'Selected %s molecule %d is not complete in the displayed supercell.\n' ...
        '  n_display_fragments       = %d\n' ...
        '  largest_fragment_fraction = %.6f'], ...
        label, row.molecule_id, ...
        row.n_display_fragments, ...
        row.largest_fragment_fraction);

    switch action
        case 'error'
            error('viz:plot_supercell_selection:IncompleteMolecule', '%s', msg);
        case 'warning'
            warning('viz:plot_supercell_selection:IncompleteMolecule', '%s', msg);
    end
end


function refUID = choose_center_complete_molecule(sys)
    T = sys.molecule_table;
    completeRows = find(T.is_complete_in_display);

    if isempty(completeRows)
        error('viz:plot_supercell_selection:NoCompleteMolecules', ...
            'No complete molecules are available in the displayed supercell.');
    end

    center = 0.5 * sum(sys.super_lattice, 1);
    com = T.com(completeRows, :);

    d2 = sum((com - center).^2, 2);
    [~, iMin] = min(d2);

    chosenRow = completeRows(iMin);
    refUID = T.molecule_id(chosenRow);
end