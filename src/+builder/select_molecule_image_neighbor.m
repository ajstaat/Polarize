function row = select_molecule_image_neighbor(Timg, varargin)
%SELECT_MOLECULE_IMAGE_NEIGHBOR
% Select one molecule-image row from a descriptor table produced by
% builder.list_molecule_images_relative_to_reference.
%
% This version supports:
%   - same-stack vs side-stack selection
%   - same/different base molecule
%   - plane-parallel vs not-plane-parallel selection using plane_angle_deg
%   - custom targets in axis_parallel / axis_perp
%
% Required:
%   Timg : table from builder.list_molecule_images_relative_to_reference
%
% Name-value options:
%   'Mode'                : 'nearest_same_stack'
%                           'nearest_side_stack'
%                           'custom'
%                           default = 'nearest_same_stack'
%
%   'BaseMolMode'         : 'any' | 'same' | 'different', default 'any'
%
%   'RequirePlaneParallel': [] | true | false
%                           []    -> ignore plane relation
%                           true  -> require plane_angle_deg <= PlaneAngleTolDeg
%                           false -> require plane_angle_deg >  PlaneAngleTolDeg
%
%   'PlaneAngleTolDeg'    : default 10
%
%   'MaxAxisPerp'         : [] or scalar
%   'MinAxisPerp'         : [] or scalar
%   'TargetAxisPerp'      : [] or scalar
%   'AxisPerpTol'         : default 0.5
%
%   'MaxAbsAxisParallel'  : [] or scalar
%   'TargetAxisParallel'  : [] or scalar
%   'AxisParallelTol'     : default 0.5
%
%   'ExcludeReference'    : default true
%
% Output:
%   row : one-row table

    p = inputParser;
    addRequired(p, 'Timg', @(x) istable(x) && ~isempty(x));

    addParameter(p, 'Mode', 'nearest_same_stack', @(x) ischar(x) || isstring(x));
    addParameter(p, 'BaseMolMode', 'any', @(x) ischar(x) || isstring(x));

    addParameter(p, 'RequirePlaneParallel', [], @(x) isempty(x) || (islogical(x) && isscalar(x)));
    addParameter(p, 'PlaneAngleTolDeg', 10, @(x) isnumeric(x) && isscalar(x) && x >= 0);

    addParameter(p, 'MaxAxisPerp', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'MinAxisPerp', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'TargetAxisPerp', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'AxisPerpTol', 0.5, @(x) isnumeric(x) && isscalar(x) && x >= 0);

    addParameter(p, 'MaxAbsAxisParallel', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'TargetAxisParallel', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'AxisParallelTol', 0.5, @(x) isnumeric(x) && isscalar(x) && x >= 0);

    addParameter(p, 'ExcludeReference', true, @(x) islogical(x) && isscalar(x));
    parse(p, Timg, varargin{:});
    opt = p.Results;

    reqVars = {'is_reference','base_mol_id','axis_parallel','axis_perp'};
    missing = setdiff(reqVars, Timg.Properties.VariableNames);
    if ~isempty(missing)
        error('Timg is missing required columns: %s', strjoin(missing, ', '));
    end

    T = Timg;

    refRows = T(T.is_reference, :);
    if height(refRows) ~= 1
        error('Timg must contain exactly one reference row.');
    end
    refBaseMol = refRows.base_mol_id(1);

    keep = true(height(T), 1);

    if opt.ExcludeReference
        keep = keep & ~T.is_reference;
    end

    switch lower(string(opt.BaseMolMode))
        case "same"
            keep = keep & (T.base_mol_id == refBaseMol);
        case "different"
            keep = keep & (T.base_mol_id ~= refBaseMol);
        case "any"
            % no-op
        otherwise
            error('BaseMolMode must be ''any'', ''same'', or ''different''.');
    end

    if ~isempty(opt.RequirePlaneParallel)
        if ~ismember('plane_angle_deg', T.Properties.VariableNames)
            error('Timg has no plane_angle_deg column. Rebuild with ComputePlaneNormals=true.');
        end

        if opt.RequirePlaneParallel
            keep = keep & (T.plane_angle_deg <= opt.PlaneAngleTolDeg);
        else
            keep = keep & (T.plane_angle_deg > opt.PlaneAngleTolDeg);
        end
    end

    if ~isempty(opt.MaxAxisPerp)
        keep = keep & (T.axis_perp <= opt.MaxAxisPerp);
    end
    if ~isempty(opt.MinAxisPerp)
        keep = keep & (T.axis_perp >= opt.MinAxisPerp);
    end
    if ~isempty(opt.MaxAbsAxisParallel)
        keep = keep & (abs(T.axis_parallel) <= opt.MaxAbsAxisParallel);
    end

    if ~isempty(opt.TargetAxisPerp)
        keep = keep & (abs(T.axis_perp - opt.TargetAxisPerp) <= opt.AxisPerpTol);
    end
    if ~isempty(opt.TargetAxisParallel)
        keep = keep & (abs(T.axis_parallel - opt.TargetAxisParallel) <= opt.AxisParallelTol);
    end

    Tcand = T(keep, :);

    if isempty(Tcand)
        error('No candidate molecule images satisfy the selector criteria.');
    end

    mode = lower(string(opt.Mode));

    switch mode
        case "nearest_same_stack"
            Tcand = Tcand(Tcand.axis_perp <= opt.AxisPerpTol, :);
            if isempty(Tcand)
                error('No same-stack candidates found within AxisPerpTol.');
            end
            [~, idx] = min(abs(Tcand.axis_parallel));
            row = Tcand(idx, :);

        case "nearest_side_stack"
            Tcand = Tcand(Tcand.axis_perp > opt.AxisPerpTol, :);
            if isempty(Tcand)
                error('No side-stack candidates found.');
            end

            d0 = min(Tcand.axis_perp);
            shell = Tcand(abs(Tcand.axis_perp - d0) <= opt.AxisPerpTol, :);

            [~, idx] = min(abs(shell.axis_parallel));
            row = shell(idx, :);

        case "custom"
            if ismember('distance', Tcand.Properties.VariableNames)
                [~, order] = sortrows([Tcand.axis_perp, abs(Tcand.axis_parallel), Tcand.distance], [1 2 3]);
            else
                [~, order] = sortrows([Tcand.axis_perp, abs(Tcand.axis_parallel)], [1 2]);
            end
            row = Tcand(order(1), :);

        otherwise
            error('Unknown Mode: %s', opt.Mode);
    end
end