function Timg = list_molecule_images_relative_to_reference(filename, varargin)
%LIST_MOLECULE_IMAGES_RELATIVE_TO_REFERENCE
% General builder for nearby molecule-image descriptors in a periodic
% molecular crystal, relative to a chosen reference molecule image.
%
% This is intended as the generalized replacement for more specialized
% stack-frame / orientation listing helpers.
%
% Required input:
%   filename : CONTCAR/POSCAR path
%
% Required name-value:
%   'RefMolID' : base molecule ID used as the reference molecule
%
% Optional name-value:
%   'RefShift'        : reference image shift [ix iy iz], default [0 0 0]
%   'BondScale'       : default 1.20
%   'SortMolecules'   : default false
%   'SearchImages'    : default [3 3 3]
%   'COMCutoff'       : default inf
%   'AxisMode'        : 'lattice' | 'plane_normal' | 'long_axis' | 'custom'
%   'AxisSpec'        : for AxisMode='lattice', use 'a'/'b'/'c' or 1/2/3
%                       for AxisMode='custom', use a 3-vector
%                       default = 'b'
%   'ComputePlaneNormals' : default true
%   'ComputeLongAxes'     : default true
%
% Output table Timg includes:
%   base_mol_id
%   ix iy iz
%   cx cy cz
%   dx dy dz
%   distance
%   axis_parallel
%   axis_perp
%   is_reference
%
% If plane normals are computed:
%   normal_x normal_y normal_z
%   plane_dot_with_ref
%   plane_angle_deg
%   orientation_sign
%
% If long axes are computed:
%   long_x long_y long_z
%   long_dot_with_ref
%   long_angle_deg
%
% Notes:
%   - This function works from whole unwrapped unit-cell molecules plus
%     lattice translations, not from wrapped sys.site_pos.
%   - It is meant to be the general descriptor table that later selectors
%     can filter however they want.

    p = inputParser;
    addRequired(p, 'filename', @(x) ischar(x) || isstring(x));

    addParameter(p, 'RefMolID', [], @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'RefShift', [0 0 0], @(x) isnumeric(x) && numel(x) == 3);
    addParameter(p, 'BondScale', 1.20, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'SortMolecules', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'SearchImages', [3 3 3], @(x) isnumeric(x) && numel(x) == 3);
    addParameter(p, 'COMCutoff', inf, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'AxisMode', 'lattice', @(x) ischar(x) || isstring(x));
    addParameter(p, 'AxisSpec', 'b');
    addParameter(p, 'ComputePlaneNormals', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ComputeLongAxes', true, @(x) islogical(x) && isscalar(x));
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

    % ------------------------------------------------------------
    % Per-base-molecule descriptors
    % ------------------------------------------------------------
    normals = nan(nBase, 3);
    longAxes = nan(nBase, 3);

    for b = 1:nBase
        X = mols{b}.cart;
        X0 = X - mean(X, 1);

        [~, ~, V] = svd(X0, 'econ');

        if opt.ComputeLongAxes
            vLong = V(:,1).';
            longAxes(b,:) = normalize_row(vLong);
        end

        if opt.ComputePlaneNormals
            vNorm = V(:,end).';
            normals(b,:) = normalize_row(vNorm);
        end
    end

    % Reference descriptors
    refMol = mols{opt.RefMolID};
    refShift = reshape(opt.RefShift, 1, 3);
    refCOM = mean((refMol.fracUnwrapped + refShift) * H, 1);

    refNormal = [nan nan nan];
    if opt.ComputePlaneNormals
        refNormal = normals(opt.RefMolID,:);
    end

    refLong = [nan nan nan];
    if opt.ComputeLongAxes
        refLong = longAxes(opt.RefMolID,:);
    end

    % ------------------------------------------------------------
    % Choose analysis axis
    % ------------------------------------------------------------
    axisVec = choose_axis_vector(H, opt.AxisMode, opt.AxisSpec, refNormal, refLong);
    axisVec = normalize_row(axisVec);

    % ------------------------------------------------------------
    % Generate nearby molecule images
    % ------------------------------------------------------------
    sx = -opt.SearchImages(1):opt.SearchImages(1);
    sy = -opt.SearchImages(2):opt.SearchImages(2);
    sz = -opt.SearchImages(3):opt.SearchImages(3);

    rows = [];

    for b = 1:nBase
        f0 = mols{b}.fracUnwrapped;

        for ix = sx
            for iy = sy
                for iz = sz
                    sh = [ix iy iz];
                    cartImg = (f0 + sh) * H;
                    comImg = mean(cartImg, 1);

                    dr = comImg - refCOM;
                    dist = norm(dr);

                    if dist <= opt.COMCutoff
                        dpar = dot(dr, axisVec);
                        dperp = norm(dr - dpar * axisVec);

                        isRef = (b == opt.RefMolID) && all(sh == refShift);

                        row = [ ...
                            b, ix, iy, iz, ...
                            comImg(1), comImg(2), comImg(3), ...
                            dr(1), dr(2), dr(3), ...
                            dist, dpar, dperp, double(isRef)];

                        if opt.ComputePlaneNormals
                            n = normals(b,:);
                            planeDot = dot(refNormal, n);
                            planeDot = max(-1, min(1, planeDot));
                            planeAng = acosd(abs(planeDot));
                            orientSign = sign(planeDot);
                            if orientSign == 0
                                orientSign = 1;
                            end
                            row = [row, n, planeDot, planeAng, orientSign];
                        end

                        if opt.ComputeLongAxes
                            u = longAxes(b,:);
                            longDot = dot(refLong, u);
                            longDot = max(-1, min(1, longDot));
                            longAng = acosd(abs(longDot));
                            row = [row, u, longDot, longAng];
                        end

                        rows = [rows; row]; %#ok<AGROW>
                    end
                end
            end
        end
    end

    % ------------------------------------------------------------
    % Table assembly
    % ------------------------------------------------------------
    varNames = { ...
        'base_mol_id','ix','iy','iz', ...
        'cx','cy','cz', ...
        'dx','dy','dz', ...
        'distance','axis_parallel','axis_perp','is_reference'};

    if opt.ComputePlaneNormals
        varNames = [varNames, ...
            {'normal_x','normal_y','normal_z','plane_dot_with_ref','plane_angle_deg','orientation_sign'}];
    end

    if opt.ComputeLongAxes
        varNames = [varNames, ...
            {'long_x','long_y','long_z','long_dot_with_ref','long_angle_deg'}];
    end

    Timg = array2table(rows, 'VariableNames', varNames);

    Timg.is_reference = logical(Timg.is_reference);

    if opt.ComputePlaneNormals
        orientation = strings(height(Timg),1);
        orientation(Timg.orientation_sign >= 0) = "parallel";
        orientation(Timg.orientation_sign < 0)  = "skew";
        Timg.orientation = orientation;
    end

    % Helpful default sorting
    Timg = sortrows(Timg, {'axis_perp','axis_parallel','distance'});
end

% -------------------------------------------------------------------------
function v = choose_axis_vector(H, axisMode, axisSpec, refNormal, refLong)
    mode = lower(string(axisMode));

    switch mode
        case "lattice"
            spec = lower(string(axisSpec));
            switch spec
                case {"a","1"}
                    v = H(1,:);
                case {"b","2"}
                    v = H(2,:);
                case {"c","3"}
                    v = H(3,:);
                otherwise
                    error('For AxisMode=''lattice'', AxisSpec must be ''a'',''b'',''c'' or 1/2/3.');
            end

        case "plane_normal"
            if any(isnan(refNormal))
                error('Plane normal not available. Set ComputePlaneNormals=true.');
            end
            v = refNormal;

        case "long_axis"
            if any(isnan(refLong))
                error('Long axis not available. Set ComputeLongAxes=true.');
            end
            v = refLong;

        case "custom"
            if ~(isnumeric(axisSpec) && numel(axisSpec) == 3)
                error('For AxisMode=''custom'', AxisSpec must be a 3-vector.');
            end
            v = reshape(axisSpec, 1, 3);

        otherwise
            error('Unknown AxisMode: %s', axisMode);
    end
end

% -------------------------------------------------------------------------
function v = normalize_row(v)
    n = norm(v);
    if n < 1e-14
        error('Cannot normalize near-zero vector.');
    end
    v = v / n;
end