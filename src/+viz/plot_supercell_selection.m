function out = plot_supercell_selection(sys, refOrResult, nbrMolID, varargin)
%PLOT_SUPERCELL_SELECTION Plot full supercell with highlighted selected pair.
%
% out = viz.plot_supercell_selection(sys, refMolID, nbrMolID)
% out = viz.plot_supercell_selection(sys, selectionResult)
% out = viz.plot_supercell_selection(..., Name, Value)
%
% This is the high-level visualization helper for Builder-selected molecule
% pairs. It shows a translucent supercell environment, highlights the
% reference and neighbor molecules, draws the supercell box, and marks the
% pair COMs and midpoint.
%
% Accepted call styles:
%
%   viz.plot_supercell_selection(sys, refMolID, nbrMolID)
%
%   result = builder.select_centered_neighbor_pair(...);
%   viz.plot_supercell_selection(sys, result)
%
% Options
%   'Axes'              axes handle, default creates one
%
%   'EnvironmentMode'   'all'|'complete'|'none', default 'all'
%                       all      : plot all environment sites
%                       complete : plot only complete molecule sites
%                       none     : plot only selected pair
%
%   'EnvironmentAlpha'  default 0.08
%   'EnvironmentSize'   default 18
%   'PairSize'          default 80
%
%   'ShowBonds'         logical, default true
%   'ShowCOM'           logical, default true
%   'ShowLabels'        logical, default true
%   'ShowMidpoint'      logical, default true
%   'DrawBox'           logical, default true
%
%   'RequireComplete'   logical, default true
%
%   'Title'             title string, default inferred
%
% Output fields include:
%   .fig
%   .ax
%   .selection_result
%   .ref_mol_id
%   .neighbor_mol_id
%   .ref_site_indices
%   .neighbor_site_indices
%   .ref_complete
%   .neighbor_complete
%   .ref_com
%   .neighbor_com
%   .pair_vector
%   .pair_distance
%   .pair_midpoint
%   .environment
%   .ref_handle
%   .neighbor_handle
%   .bond_handle
%   .com_handle
%   .midpoint_handle
%   .pair_line
%   .label_handle

if nargin < 2
    error('viz:plot_supercell_selection:NotEnoughInputs', ...
        'Expected sys plus either a selection result or ref/nbr molecule IDs.');
end

if nargin < 3
    nbrMolID = [];
end

[refMolID, nbrMolID, selectionResult] = local_parse_selection_input(refOrResult, nbrMolID);

p = inputParser;
addRequired(p, 'sys', @isstruct);
addRequired(p, 'refOrResult');
addOptional(p, 'nbrMolID', nbrMolID);
addParameter(p, 'Axes', [], @(x) isempty(x) || isgraphics(x, 'axes'));
addParameter(p, 'EnvironmentMode', 'all', @(x) ischar(x) || isstring(x));
addParameter(p, 'EnvironmentAlpha', 0.08, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'EnvironmentSize', 18, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'PairSize', 80, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'ShowBonds', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ShowCOM', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ShowLabels', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ShowMidpoint', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'DrawBox', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'RequireComplete', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));
parse(p, sys, refOrResult, nbrMolID, varargin{:});

opt = p.Results;

validate_sys(sys);

if opt.RequireComplete
    local_assert_complete(sys, [refMolID nbrMolID]);
end

[refComplete, nbrComplete] = local_pair_completeness(sys, refMolID, nbrMolID);

idxRef = builder.site_indices_for_molecule(sys, refMolID);
idxNbr = builder.site_indices_for_molecule(sys, nbrMolID);

if isempty(idxRef)
    error('viz:plot_supercell_selection:UnknownReferenceMolecule', ...
        'Reference molecule ID %d was not found.', refMolID);
end

if isempty(idxNbr)
    error('viz:plot_supercell_selection:UnknownNeighborMolecule', ...
        'Neighbor molecule ID %d was not found.', nbrMolID);
end

Xref = sys.site_pos(idxRef, :);
Xnbr = sys.site_pos(idxNbr, :);

refCOM = mean(Xref, 1);
nbrCOM = mean(Xnbr, 1);

pairVector = nbrCOM - refCOM;
pairDistance = norm(pairVector);
pairMidpoint = 0.5 * (refCOM + nbrCOM);

ax = opt.Axes;

if isempty(ax)
    fig = figure;
    ax = axes(fig);
else
    fig = ancestor(ax, 'figure');
end

holdState = ishold(ax);
hold(ax, 'on');

% -------------------------------------------------------------------------
% Environment/background
% -------------------------------------------------------------------------

environmentMode = lower(char(string(opt.EnvironmentMode)));

env = struct();
env.fig = fig;
env.ax = ax;
env.site_handle = gobjects(0);
env.box_handle = gobjects(0);
env.com_handle = gobjects(0);
env.bond_handle = gobjects(0);
env.site_indices = [];

switch environmentMode
    case 'all'
        env = viz.plot_system(sys, ...
            'Axes', ax, ...
            'ColorBy', 'complete', ...
            'MarkerSize', opt.EnvironmentSize, ...
            'Alpha', opt.EnvironmentAlpha, ...
            'DrawBox', opt.DrawBox, ...
            'ShowCOM', false, ...
            'ShowBonds', false);

    case 'complete'
        completeIDs = builder.complete_molecule_ids(sys);

        env = viz.plot_system(sys, ...
            'Axes', ax, ...
            'MoleculeIDs', completeIDs, ...
            'ColorBy', 'molecule', ...
            'MarkerSize', opt.EnvironmentSize, ...
            'Alpha', opt.EnvironmentAlpha, ...
            'DrawBox', opt.DrawBox, ...
            'ShowCOM', false, ...
            'ShowBonds', false);

    case 'none'
        if opt.DrawBox
            if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
                env.box_handle = viz.draw_cell_box(ax, sys.super_lattice, ...
                    'VectorsAs', 'rows', ...
                    'Color', [0 0 0], ...
                    'LineWidth', 1.0);
            elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
                env.box_handle = viz.draw_cell_box(ax, sys.lattice, ...
                    'VectorsAs', 'rows', ...
                    'Color', [0 0 0], ...
                    'LineWidth', 1.0);
            end
        end

    otherwise
        error('viz:plot_supercell_selection:BadEnvironmentMode', ...
            'EnvironmentMode must be ''all'', ''complete'', or ''none''.');
end

% -------------------------------------------------------------------------
% Selected pair overlay
% -------------------------------------------------------------------------

refHandle = local_scatter_selected(ax, Xref, opt.PairSize, 'ref');
nbrHandle = local_scatter_selected(ax, Xnbr, opt.PairSize, 'nbr');

bondHandle = gobjects(0);

if opt.ShowBonds
    bondHandle = [
        local_plot_molecule_bonds(ax, sys, idxRef)
        local_plot_molecule_bonds(ax, sys, idxNbr)
    ];
end

comHandle = gobjects(0);
pairLine = gobjects(0);
midpointHandle = gobjects(0);

if opt.ShowCOM
    comHandle = scatter3(ax, ...
        [refCOM(1); nbrCOM(1)], ...
        [refCOM(2); nbrCOM(2)], ...
        [refCOM(3); nbrCOM(3)], ...
        120, 'kx', 'LineWidth', 2);

    pairLine = plot3(ax, ...
        [refCOM(1) nbrCOM(1)], ...
        [refCOM(2) nbrCOM(2)], ...
        [refCOM(3) nbrCOM(3)], ...
        'k--', 'LineWidth', 1.25);
end

if opt.ShowMidpoint
    midpointHandle = scatter3(ax, ...
        pairMidpoint(1), pairMidpoint(2), pairMidpoint(3), ...
        90, 'ko', 'LineWidth', 1.5);
end

labelHandle = gobjects(0);

if opt.ShowLabels
    labelHandle = gobjects(2, 1);

    labelHandle(1) = text(ax, refCOM(1), refCOM(2), refCOM(3), ...
        sprintf('  ref %d', refMolID), ...
        'FontWeight', 'bold', ...
        'Color', [0 0 0], ...
        'BackgroundColor', [1 1 1], ...
        'Margin', 1);

    labelHandle(2) = text(ax, nbrCOM(1), nbrCOM(2), nbrCOM(3), ...
        sprintf('  nbr %d', nbrMolID), ...
        'FontWeight', 'bold', ...
        'Color', [0 0 0], ...
        'BackgroundColor', [1 1 1], ...
        'Margin', 1);
end

axis(ax, 'equal');
grid(ax, 'on');
xlabel(ax, 'x / bohr');
ylabel(ax, 'y / bohr');
zlabel(ax, 'z / bohr');
view(ax, 3);

ttl = char(string(opt.Title));

if isempty(ttl)
    ttl = local_default_title(selectionResult, refMolID, nbrMolID);
end

title(ax, ttl);

if ~holdState
    hold(ax, 'off');
end

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------

out = struct();

out.fig = fig;
out.ax = ax;

out.selection_result = selectionResult;

out.ref_mol_id = refMolID;
out.neighbor_mol_id = nbrMolID;

out.ref_site_indices = idxRef(:);
out.neighbor_site_indices = idxNbr(:);

out.ref_complete = refComplete;
out.neighbor_complete = nbrComplete;

out.ref_com = refCOM;
out.neighbor_com = nbrCOM;

out.pair_vector = pairVector;
out.pair_distance = pairDistance;
out.pair_midpoint = pairMidpoint;

out.environment_mode = environmentMode;
out.environment = env;

out.ref_handle = refHandle;
out.neighbor_handle = nbrHandle;

out.bond_handle = bondHandle;
out.com_handle = comHandle;
out.midpoint_handle = midpointHandle;
out.pair_line = pairLine;
out.label_handle = labelHandle;

end

% =========================================================================
% Input parsing
% =========================================================================

function [refMolID, nbrMolID, selectionResult] = local_parse_selection_input(refOrResult, nbrMolID)

selectionResult = [];

if isstruct(refOrResult)
    selectionResult = refOrResult;

    if isfield(selectionResult, 'reference_mol_id')
        refMolID = selectionResult.reference_mol_id;
    elseif isfield(selectionResult, 'ref_mol_id')
        refMolID = selectionResult.ref_mol_id;
    else
        error('viz:plot_supercell_selection:BadSelectionResult', ...
            'Selection result must contain reference_mol_id or ref_mol_id.');
    end

    if isfield(selectionResult, 'neighbor_mol_id')
        nbrMolID = selectionResult.neighbor_mol_id;
    elseif isfield(selectionResult, 'nbr_mol_id')
        nbrMolID = selectionResult.nbr_mol_id;
    else
        error('viz:plot_supercell_selection:BadSelectionResult', ...
            'Selection result must contain neighbor_mol_id or nbr_mol_id.');
    end
else
    refMolID = refOrResult;

    if isempty(nbrMolID)
        error('viz:plot_supercell_selection:MissingNeighborMolID', ...
            'Neighbor molecule ID is required when not passing a selection result struct.');
    end
end

if ~(isnumeric(refMolID) && isscalar(refMolID) && isfinite(refMolID))
    error('viz:plot_supercell_selection:BadReferenceMolID', ...
        'Reference molecule ID must be a finite numeric scalar.');
end

if ~(isnumeric(nbrMolID) && isscalar(nbrMolID) && isfinite(nbrMolID))
    error('viz:plot_supercell_selection:BadNeighborMolID', ...
        'Neighbor molecule ID must be a finite numeric scalar.');
end

end

% =========================================================================
% Validation / completeness
% =========================================================================

function validate_sys(sys)

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos) || size(sys.site_pos, 2) ~= 3
    error('viz:plot_supercell_selection:BadSitePos', ...
        'sys.site_pos must be N x 3.');
end

n = size(sys.site_pos, 1);

if ~isfield(sys, 'site_mol_id') || isempty(sys.site_mol_id) || numel(sys.site_mol_id) ~= n
    error('viz:plot_supercell_selection:BadSiteMolID', ...
        'sys.site_mol_id must have one entry per site.');
end

if ~isfield(sys, 'site_type') || isempty(sys.site_type) || numel(sys.site_type) ~= n
    error('viz:plot_supercell_selection:BadSiteType', ...
        'sys.site_type must have one entry per site.');
end

if ~isfield(sys, 'molecule_table') || isempty(sys.molecule_table)
    error('viz:plot_supercell_selection:MissingMoleculeTable', ...
        'sys.molecule_table is required.');
end

T = sys.molecule_table;

if ~isfield(T, 'molecule_id') || isempty(T.molecule_id) || ...
        ~isfield(T, 'is_complete_in_display') || isempty(T.is_complete_in_display)
    error('viz:plot_supercell_selection:BadMoleculeTable', ...
        'sys.molecule_table must contain molecule_id and is_complete_in_display.');
end

end

function local_assert_complete(sys, molIDs)

T = sys.molecule_table;

for k = 1:numel(molIDs)
    row = find(T.molecule_id == molIDs(k), 1, 'first');

    if isempty(row)
        error('viz:plot_supercell_selection:UnknownMolecule', ...
            'Molecule ID %d was not found.', molIDs(k));
    end

    if ~T.is_complete_in_display(row)
        error('viz:plot_supercell_selection:IncompleteMolecule', ...
            'Molecule ID %d is not complete in the displayed supercell.', molIDs(k));
    end
end

end

function [refComplete, nbrComplete] = local_pair_completeness(sys, refMolID, nbrMolID)

T = sys.molecule_table;

refRow = find(T.molecule_id == refMolID, 1, 'first');
nbrRow = find(T.molecule_id == nbrMolID, 1, 'first');

if isempty(refRow)
    error('viz:plot_supercell_selection:UnknownReferenceMolecule', ...
        'Reference molecule ID %d was not found.', refMolID);
end

if isempty(nbrRow)
    error('viz:plot_supercell_selection:UnknownNeighborMolecule', ...
        'Neighbor molecule ID %d was not found.', nbrMolID);
end

refComplete = logical(T.is_complete_in_display(refRow));
nbrComplete = logical(T.is_complete_in_display(nbrRow));

end

% =========================================================================
% Plot helpers
% =========================================================================

function h = local_scatter_selected(ax, X, markerSize, role)

switch role
    case 'ref'
        marker = 'o';

    case 'nbr'
        marker = '^';

    otherwise
        marker = 'o';
end

h = scatter3(ax, X(:,1), X(:,2), X(:,3), ...
    markerSize, marker, ...
    'filled', ...
    'MarkerEdgeColor', [0 0 0], ...
    'LineWidth', 0.75);

end

function h = local_plot_molecule_bonds(ax, sys, idx)

X = sys.site_pos(idx, :);
types = local_to_cell_column(sys.site_type(idx));
toAng = local_length_to_angstrom_factor(sys);

n = numel(idx);
h = gobjects(0);
count = 0;

for i = 1:(n - 1)
    for j = (i + 1):n
        dAng = norm(X(j,:) - X(i,:)) * toAng;

        if io.bond_graph_tools.is_bonded(types{i}, types{j}, dAng, 1.20)
            count = count + 1;
            pts = [X(i,:); X(j,:)];
            h(count,1) = plot3(ax, pts(:,1), pts(:,2), pts(:,3), ...
                '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.1); %#ok<AGROW>
        end
    end
end

end

% =========================================================================
% Title / units / text helpers
% =========================================================================

function ttl = local_default_title(selectionResult, refMolID, nbrMolID)

if isempty(selectionResult)
    ttl = sprintf('Selected pair: ref %d, neighbor %d', refMolID, nbrMolID);
    return;
end

relation = '';

if isfield(selectionResult, 'relation')
    relation = char(string(selectionResult.relation));
elseif isfield(selectionResult, 'neighbor_type')
    relation = char(string(selectionResult.neighbor_type));
end

shellText = '';

if isfield(selectionResult, 'shell')
    shellText = sprintf(' shell %d', selectionResult.shell);
elseif isfield(selectionResult, 'requested_shell')
    shellText = sprintf(' shell %d', selectionResult.requested_shell);
end

if isempty(relation)
    ttl = sprintf('Selected pair: ref %d, neighbor %d', refMolID, nbrMolID);
else
    ttl = sprintf('%s%s: ref %d, neighbor %d', ...
        relation, shellText, refMolID, nbrMolID);
end

end

function factor = local_length_to_angstrom_factor(sys)

BOHR2ANG = 1 / 1.8897259886;

if isfield(sys, 'units') && isfield(sys.units, 'length')
    unit = lower(char(string(sys.units.length)));
else
    unit = 'bohr';
end

switch unit
    case {'angstrom', 'angstroms', 'ang', 'a'}
        factor = 1.0;

    case {'bohr', 'bohrs', 'a0', 'au_length'}
        factor = BOHR2ANG;

    otherwise
        error('viz:plot_supercell_selection:UnsupportedLengthUnit', ...
            'Unsupported length unit: %s', unit);
end

end

function c = local_to_cell_column(x)

if isstring(x)
    c = cellstr(x(:));
elseif iscellstr(x)
    c = x(:);
elseif iscell(x)
    c = x(:);
else
    error('viz:plot_supercell_selection:BadTextMetadata', ...
        'Text metadata must be string/cellstr/cell.');
end

end