function info = molecule_display_status(sys, uniqueMolID, varargin)
%MOLECULE_DISPLAY_STATUS Check whether one supercell molecule is complete/
% contiguous in the displayed supercell box.
%
% info = builder.molecule_display_status(sys, uniqueMolID)
% info = builder.molecule_display_status(sys, uniqueMolID, 'BondScale', 1.20)
%
% Inputs
%   sys          working system struct with at least:
%                  .site_pos
%                  .site_type
%                  .site_mol_id
%   uniqueMolID  scalar supercell molecule ID
%
% Optional name-value inputs
%   'BondScale'  default 1.20
%
% Output
%   info struct with fields:
%       .unique_mol_id
%       .site_indices
%       .n_sites
%       .fragment_ids
%       .fragment_sizes
%       .n_fragments
%       .is_complete_in_display
%       .largest_fragment_fraction
%       .com

    p = inputParser;
    addRequired(p, 'sys', @isstruct);
    addRequired(p, 'uniqueMolID', @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'BondScale', 1.20, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, sys, uniqueMolID, varargin{:});
    opt = p.Results;

    validate_sys(sys);

    idx = builder.site_indices_for_molecule(sys, uniqueMolID);
    if isempty(idx)
        error('builder:molecule_display_status:UnknownMolID', ...
            'uniqueMolID %d was not found in sys.site_mol_id.', uniqueMolID);
    end

    X = sys.site_pos(idx, :);
    species = sys.site_type(idx);

    n = size(X, 1);
    if n == 1
        fragID = 1;
    else
        A = false(n, n);

        for i = 1:n-1
            si = species{i};

            for j = i+1:n
                sj = species{j};
                d = norm(X(i,:) - X(j,:));

                if io.bond_graph_tools.is_bonded(si, sj, d, opt.BondScale)
                    A(i,j) = true;
                    A(j,i) = true;
                end
            end
        end

        fragID = conncomp(graph(A)).';
    end

    fragIDs = unique(fragID, 'stable');
    nFrag = numel(fragIDs);
    fragSizes = zeros(nFrag, 1);

    for k = 1:nFrag
        fragSizes(k) = nnz(fragID == fragIDs(k));
    end

    info = struct();
    info.unique_mol_id = uniqueMolID;
    info.site_indices = idx;
    info.n_sites = n;
    info.fragment_ids = fragIDs;
    info.fragment_sizes = fragSizes;
    info.n_fragments = nFrag;
    info.is_complete_in_display = (nFrag == 1);
    info.largest_fragment_fraction = max(fragSizes) / n;
    info.com = mean(X, 1);
end

function validate_sys(sys)
    required = {'site_pos', 'site_type', 'site_mol_id'};
    for k = 1:numel(required)
        name = required{k};
        if ~isfield(sys, name) || isempty(sys.(name))
            error('builder:molecule_display_status:MissingField', ...
                'sys.%s is required and missing/empty.', name);
        end
    end

    if size(sys.site_pos, 2) ~= 3
        error('builder:molecule_display_status:BadSitePos', ...
            'sys.site_pos must be N x 3.');
    end

    if numel(sys.site_type) ~= size(sys.site_pos, 1)
        error('builder:molecule_display_status:BadSiteType', ...
            'sys.site_type must have one entry per site.');
    end
end