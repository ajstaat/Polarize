function out = plot_supercell_selection(sys, refOrResult, nbrMolID, varargin)
%PLOT_SUPERCELL_SELECTION Plot a selected molecule pair in a supercell.
%
% Style target:
%   - light grey environment
%   - bold highlighted pair
%   - no bonds
%   - optional COM markers, labels, and pair line
%   - explicit cell box
%
% Usage
%   out = viz.plot_supercell_selection(sys, refMolID, nbrMolID)
%   out = viz.plot_supercell_selection(sys, selectionResult)
%   out = viz.plot_supercell_selection(..., Name, Value)
%
% Accepted selection-result fields:
%   .reference_mol_id or .ref_mol_id
%   .neighbor_mol_id  or .nbr_mol_id
%
% Options
%   'Axes'         axes handle, default creates one
%   'Title'        plot title, default inferred
%   'DrawBox'      logical, default true
%   'ShowCOM'      logical, default true
%   'ShowLabels'   logical, default true
%   'ShowPairLine' logical, default true
%
% Output
%   out.fig
%   out.ax
%   out.ref_mol_id
%   out.neighbor_mol_id
%   out.ref_site_indices
%   out.neighbor_site_indices
%   out.ref_com
%   out.neighbor_com
%   out.pair_vector
%   out.pair_midpoint
%   out.pair_distance
%   out.environment_handle
%   out.ref_handle
%   out.neighbor_handle
%   out.com_handle
%   out.label_handle
%   out.pair_line
%   out.box_handle

if nargin < 2
    error('viz:plot_supercell_selection:NotEnoughInputs', ...
        'Expected sys plus either a selection struct or ref/nbr molecule IDs.');
end

if nargin < 3
    nbrMolID = [];
end

% Support:
%   viz.plot_supercell_selection(sys, selectionStruct, 'Axes', ax, ...)
if isstruct(refOrResult) && (ischar(nbrMolID) || isstring(nbrMolID))
    varargin = [{nbrMolID}, varargin];
    nbrMolID = [];
end

[refMolID, nbrMolID] = local_parse_pair(refOrResult, nbrMolID);

p = inputParser;
addRequired(p, 'sys', @isstruct);
addRequired(p, 'refOrResult');
addOptional(p, 'nbrMolID', nbrMolID);
addParameter(p, 'Axes', [], @(x) isempty(x) || isgraphics(x, 'axes'));
addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'DrawBox', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ShowCOM', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ShowLabels', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ShowPairLine', true, @(x) islogical(x) && isscalar(x));
parse(p, sys, refOrResult, nbrMolID, varargin{:});

opt = p.Results;

local_validate_sys(sys);

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

Xall = sys.site_pos;
Xref = sys.site_pos(idxRef, :);
Xnbr = sys.site_pos(idxNbr, :);

refCOM = mean(Xref, 1);
nbrCOM = mean(Xnbr, 1);

pairVector = nbrCOM - refCOM;
pairMidpoint = 0.5 * (refCOM + nbrCOM);
pairDistance = norm(pairVector);

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
% Environment: large pale-grey markers
% -------------------------------------------------------------------------

environmentHandle = scatter3(ax, ...
    Xall(:,1), Xall(:,2), Xall(:,3), ...
    36, ...
    'o', ...
    'filled', ...
    'MarkerFaceColor', [0.75 0.75 0.75], ...
    'MarkerEdgeColor', [0.75 0.75 0.75]);

try
    environmentHandle.MarkerFaceAlpha = 0.22;
    environmentHandle.MarkerEdgeAlpha = 0.22;
catch
end

% -------------------------------------------------------------------------
% Highlighted pair: bold markers
% -------------------------------------------------------------------------

refHandle = scatter3(ax, ...
    Xref(:,1), Xref(:,2), Xref(:,3), ...
    80, ...
    'o', ...
    'filled', ...
    'MarkerFaceColor', [0.85 0.20 0.20], ...
    'MarkerEdgeColor', [0 0 0], ...
    'LineWidth', 0.75);

nbrHandle = scatter3(ax, ...
    Xnbr(:,1), Xnbr(:,2), Xnbr(:,3), ...
    80, ...
    'o', ...
    'filled', ...
    'MarkerFaceColor', [0.20 0.35 0.85], ...
    'MarkerEdgeColor', [0 0 0], ...
    'LineWidth', 0.75);

% -------------------------------------------------------------------------
% Optional cell box
% -------------------------------------------------------------------------

boxHandle = gobjects(0);

if opt.DrawBox
    if isfield(sys, 'super_lattice') && ~isempty(sys.super_lattice)
        boxHandle = viz.draw_cell_box(ax, sys.super_lattice, ...
            'VectorsAs', 'rows', ...
            'Color', [0 0 0], ...
            'LineWidth', 1.0);
    elseif isfield(sys, 'lattice') && ~isempty(sys.lattice)
        boxHandle = viz.draw_cell_box(ax, sys.lattice, ...
            'VectorsAs', 'rows', ...
            'Color', [0 0 0], ...
            'LineWidth', 1.0);
    end
end

% -------------------------------------------------------------------------
% Optional COMs / pair line / labels
% -------------------------------------------------------------------------

comHandle = gobjects(0);
pairLine = gobjects(0);
labelHandle = gobjects(0);

if opt.ShowCOM
    comHandle = scatter3(ax, ...
        [refCOM(1); nbrCOM(1)], ...
        [refCOM(2); nbrCOM(2)], ...
        [refCOM(3); nbrCOM(3)], ...
        90, ...
        'kx', ...
        'LineWidth', 1.75);
end

if opt.ShowPairLine
    pairLine = plot3(ax, ...
        [refCOM(1) nbrCOM(1)], ...
        [refCOM(2) nbrCOM(2)], ...
        [refCOM(3) nbrCOM(3)], ...
        'k--', ...
        'LineWidth', 1.1);
end

if opt.ShowLabels
    labelHandle = gobjects(2,1);

    labelHandle(1) = text(ax, ...
        refCOM(1), refCOM(2), refCOM(3), ...
        sprintf('  ref %d', refMolID), ...
        'FontWeight', 'bold', ...
        'Color', [0 0 0]);

    labelHandle(2) = text(ax, ...
        nbrCOM(1), nbrCOM(2), nbrCOM(3), ...
        sprintf('  nbr %d', nbrMolID), ...
        'FontWeight', 'bold', ...
        'Color', [0 0 0]);
end

axis(ax, 'equal');
grid(ax, 'on');
xlabel(ax, 'x / bohr');
ylabel(ax, 'y / bohr');
zlabel(ax, 'z / bohr');
view(ax, 3);

ttl = char(string(opt.Title));
if isempty(ttl)
    ttl = sprintf('Selected pair: ref %d, neighbor %d', refMolID, nbrMolID);
end
title(ax, ttl);

if ~holdState
    hold(ax, 'off');
end

out = struct();
out.fig = fig;
out.ax = ax;

out.ref_mol_id = refMolID;
out.neighbor_mol_id = nbrMolID;

out.ref_site_indices = idxRef(:);
out.neighbor_site_indices = idxNbr(:);

out.ref_com = refCOM;
out.neighbor_com = nbrCOM;
out.pair_vector = pairVector;
out.pair_midpoint = pairMidpoint;
out.pair_distance = pairDistance;

out.environment_handle = environmentHandle;
out.ref_handle = refHandle;
out.neighbor_handle = nbrHandle;
out.com_handle = comHandle;
out.label_handle = labelHandle;
out.pair_line = pairLine;
out.box_handle = boxHandle;

end

% =========================================================================
% Helpers
% =========================================================================

function [refMolID, nbrMolID] = local_parse_pair(refOrResult, nbrMolID)

if isstruct(refOrResult)
    S = refOrResult;

    if isfield(S, 'reference_mol_id')
        refMolID = S.reference_mol_id;
    elseif isfield(S, 'ref_mol_id')
        refMolID = S.ref_mol_id;
    else
        error('viz:plot_supercell_selection:BadSelectionStruct', ...
            'Selection struct must contain reference_mol_id or ref_mol_id.');
    end

    if isfield(S, 'neighbor_mol_id')
        nbrMolID = S.neighbor_mol_id;
    elseif isfield(S, 'nbr_mol_id')
        nbrMolID = S.nbr_mol_id;
    else
        error('viz:plot_supercell_selection:BadSelectionStruct', ...
            'Selection struct must contain neighbor_mol_id or nbr_mol_id.');
    end

else
    refMolID = refOrResult;

    if isempty(nbrMolID)
        error('viz:plot_supercell_selection:MissingNeighborMolID', ...
            'Neighbor molecule ID is required.');
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

function local_validate_sys(sys)

if ~isfield(sys, 'site_pos') || isempty(sys.site_pos) || size(sys.site_pos,2) ~= 3
    error('viz:plot_supercell_selection:BadSitePos', ...
        'sys.site_pos must be N x 3.');
end

n = size(sys.site_pos, 1);

if ~isfield(sys, 'site_mol_id') || isempty(sys.site_mol_id) || numel(sys.site_mol_id) ~= n
    error('viz:plot_supercell_selection:BadSiteMolID', ...
        'sys.site_mol_id must have one entry per site.');
end

end