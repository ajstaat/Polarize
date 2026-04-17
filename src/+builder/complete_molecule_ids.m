function molIDs = complete_molecule_ids(sys, varargin)
%COMPLETE_MOLECULE_IDS Return molecule IDs that are complete in display.
%
% molIDs = builder.complete_molecule_ids(sys)
% molIDs = builder.complete_molecule_ids(sys, 'Verbose', true)
%
% Output
%   molIDs   column vector of supercell molecule IDs whose
%            molecule_table.is_complete_in_display is true

    p = inputParser;
    addRequired(p, 'sys', @isstruct);
    addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
    parse(p, sys, varargin{:});
    opt = p.Results;

    validate_sys(sys);

    T = sys.molecule_table;
    mask = T.is_complete_in_display(:);

    molIDs = T.molecule_id(mask);
    molIDs = molIDs(:);

    if opt.Verbose
        nTotal = numel(T.molecule_id);
        nComplete = numel(molIDs);

        fprintf('Complete molecule summary:\n');
        fprintf('  total molecules found       = %d\n', nTotal);
        fprintf('  complete in displayed box   = %d\n', nComplete);

        if nComplete > 0
            fprintf('  complete molecule IDs       = %s\n', mat2str(molIDs(:).'));
        else
            fprintf('  complete molecule IDs       = []\n');
        end
    end
end


function validate_sys(sys)
    if ~isstruct(sys)
        error('builder:complete_molecule_ids:BadInput', ...
            'sys must be a struct.');
    end

    if ~isfield(sys, 'molecule_table') || isempty(sys.molecule_table)
        error('builder:complete_molecule_ids:MissingMoleculeTable', ...
            'sys.molecule_table is required and missing/empty.');
    end

    T = sys.molecule_table;
    required = {'molecule_id', 'is_complete_in_display'};

    for k = 1:numel(required)
        name = required{k};
        if ~isfield(T, name) || isempty(T.(name))
            error('builder:complete_molecule_ids:BadMoleculeTable', ...
                'sys.molecule_table.%s is required and missing/empty.', name);
        end
    end
end