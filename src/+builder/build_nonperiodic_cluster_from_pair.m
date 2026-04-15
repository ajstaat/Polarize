function out = build_nonperiodic_cluster_from_pair(filename, refMolID, neighborRow, varargin)
%BUILD_NONPERIODIC_CLUSTER_FROM_PAIR
% Build a finite nonperiodic molecular cluster around a selected pair:
%   - one reference base molecule at RefShift
%   - one selected neighbor image from NeighborRow
%
% The cluster consists of whole molecule images whose COM lies within
% ClusterCutoff of the pair midpoint.
%
% Output:
%   out.crystal      : Polarize-style crystal struct for nonperiodic use
%   out.image_table  : table of included molecule images
%   out.focus_pair   : struct with ref / neighbor cluster mol IDs

    p = inputParser;
    addRequired(p, 'filename', @(x) ischar(x) || isstring(x));
    addRequired(p, 'refMolID', @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addRequired(p, 'neighborRow', @(x) istable(x) && height(x) == 1);

    addParameter(p, 'RefShift', [0 0 0], @(x) isnumeric(x) && numel(x) == 3);
    addParameter(p, 'BondScale', 1.20, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'SortMolecules', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'SearchImages', [4 4 4], @(x) isnumeric(x) && numel(x) == 3);
    addParameter(p, 'ClusterCutoff', 20.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'CenterMode', 'pair_midpoint', @(x) ischar(x) || isstring(x));
    addParameter(p, 'BoxPadding', 5.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, filename, refMolID, neighborRow, varargin{:});
    opt = p.Results;

    S = io.read_vasp_structure(filename);
    H = S.lattice;

    mols = io.unwrap_all_contcar_molecules(filename, ...
        'BondScale', opt.BondScale, ...
        'SortMolecules', opt.SortMolecules);

    refShift = reshape(opt.RefShift, 1, 3);
    refMol = mols{refMolID};
    refCart = (refMol.fracUnwrapped + refShift) * H;
    refCOM = mean(refCart, 1);

    neighBase = neighborRow.base_mol_id(1);
    neighShift = [neighborRow.ix(1), neighborRow.iy(1), neighborRow.iz(1)];
    neighMol = mols{neighBase};
    neighCart = (neighMol.fracUnwrapped + neighShift) * H;
    neighCOM = mean(neighCart, 1);

    switch lower(string(opt.CenterMode))
        case "pair_midpoint"
            center = 0.5 * (refCOM + neighCOM);
        case "reference_com"
            center = refCOM;
        case "none"
            center = [0 0 0];
        otherwise
            error('Unknown CenterMode: %s', opt.CenterMode);
    end

    sx = -opt.SearchImages(1):opt.SearchImages(1);
    sy = -opt.SearchImages(2):opt.SearchImages(2);
    sz = -opt.SearchImages(3):opt.SearchImages(3);

    rows = [];
    molCounter = 0;

    cartAll = [];
    typeAll = {};
    labelAll = {};
    molIDAll = [];

    for b = 1:numel(mols)
        Mol = mols{b};
        for ix = sx
            for iy = sy
                for iz = sz
                    sh = [ix iy iz];
                    cartImg = (Mol.fracUnwrapped + sh) * H;
                    comImg = mean(cartImg, 1);

                    if norm(comImg - center) <= opt.ClusterCutoff
                        molCounter = molCounter + 1;

                        rows = [rows; ...
                            molCounter, b, ix, iy, iz, ...
                            comImg(1), comImg(2), comImg(3)]; %#ok<AGROW>

                        n = size(cartImg, 1);
                        cartAll = [cartAll; cartImg]; %#ok<AGROW>
                        typeAll = [typeAll; S.species(Mol.indices)]; %#ok<AGROW>
                        % Build molecule-local unique labels like C1, C2, H1, ...
                        speciesThis = S.species(Mol.indices);
                        labelsThis = make_local_labels(speciesThis);
                        labelAll = [labelAll; labelsThis(:)]; %#ok<AGROW>
                        molIDAll = [molIDAll; molCounter * ones(n,1)]; %#ok<AGROW>
                    end
                end
            end
        end
    end

    % Recenter Cartesian coordinates for numerical convenience
    if ~strcmpi(opt.CenterMode, 'none')
        cartAll = cartAll - center;
        rows(:,6:8) = rows(:,6:8) - center;
    end

    Timg = array2table(rows, 'VariableNames', ...
        {'unique_mol_id','base_mol_id','ix','iy','iz','cx','cy','cz'});

    % Locate focus pair inside cluster mol IDs
    refCluster = Timg( ...
        Timg.base_mol_id == refMolID & ...
        Timg.ix == refShift(1) & ...
        Timg.iy == refShift(2) & ...
        Timg.iz == refShift(3), :);

    neighCluster = Timg( ...
        Timg.base_mol_id == neighBase & ...
        Timg.ix == neighShift(1) & ...
        Timg.iy == neighShift(2) & ...
        Timg.iz == neighShift(3), :);

    if height(refCluster) ~= 1 || height(neighCluster) ~= 1
        error('Could not uniquely locate reference or neighbor inside cluster.');
    end

    % Build a simple orthorhombic bounding box lattice for nonperiodic use
    xyzMin = min(cartAll, [], 1) - opt.BoxPadding;
    xyzMax = max(cartAll, [], 1) + opt.BoxPadding;
    boxSize = xyzMax - xyzMin;

    cartBox = cartAll - xyzMin;   % shift into positive box
    lattice = diag(boxSize);
    fracBox = cartBox / lattice;  % since lattice is diagonal rows, mrdivide is fine

    crystal = struct();
    crystal.cellpar = [];
    crystal.lattice = lattice;
    crystal.frac_coords = fracBox;
    crystal.cart_coords = cartBox;
    crystal.mol_id = molIDAll;
    crystal.site_type = typeAll(:);
    crystal.site_label = labelAll(:);

    out = struct();
    out.crystal = crystal;
    out.image_table = Timg;
    out.focus_pair = struct( ...
        'ref_unique_mol_id', refCluster.unique_mol_id, ...
        'neighbor_unique_mol_id', neighCluster.unique_mol_id, ...
        'ref_base_mol_id', refMolID, ...
        'neighbor_base_mol_id', neighBase, ...
        'ref_shift', refShift, ...
        'neighbor_shift', neighShift, ...
        'box_min_cart', xyzMin, ...
        'box_max_cart', xyzMax);

    
end
function labels = make_local_labels(species)
    %MAKE_LOCAL_LABELS Generate labels like C1, C2, H1... within one molecule.

    species = species(:);
    labels = cell(size(species));

    counts = containers.Map('KeyType','char', 'ValueType','double');

    for k = 1:numel(species)
        sym = char(species{k});

        if isKey(counts, sym)
            counts(sym) = counts(sym) + 1;
        else
            counts(sym) = 1;
        end

        labels{k} = sprintf('%s%d', sym, counts(sym));
    end
end