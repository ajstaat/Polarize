function sys = apply_molecule_charges(sys, molIDs, varargin)
%APPLY_MOLECULE_CHARGES Assign charges to one or more selected molecules.
%
% sys = builder.apply_molecule_charges(sys, molIDs, ...)
%
% Inputs
%   sys     working system struct from builder.make_crystal_system
%   molIDs  scalar or vector of supercell molecule IDs
%
% Charge modes
%   1) Uniform total charge:
%        'Mode', 'uniform'
%        'TotalCharges', [q1 q2 ...]
%
%   2) Charges from file in template order:
%        'Mode', 'file'
%        'ChargeFiles', {'mol1.txt', 'mol2.txt', ...}
%        'TemplateCoords', {X1, X2, ...}
%        'TemplateSpecies', {S1, S2, ...}
%
%      In file mode, charges are assumed to be listed in TEMPLATE atom order.
%      The helper maps template atoms onto the selected supercell molecule
%      using builder.match_molecule_atoms_by_frame.
%
% Optional name-value inputs
%   'Mode'                          'uniform' | 'file', default 'uniform'
%   'TotalCharges'                  numeric vector, required for uniform mode
%   'ChargeFiles'                   cellstr, required for file mode
%   'TemplateCoords'                cell array of N x 3 arrays, required for file mode
%   'TemplateSpecies'               cell array of species cell arrays, required for file mode
%   'ReferenceAxis'                 1x3 vector, default sys.super_lattice(2,:)
%   'DistanceTol'                   scalar, default 1e-2
%   'SetActive'                     logical, default true
%   'DisablePolarizabilityOnCharged' logical, default true
%   'ZeroExistingCharges'           logical, default true
%   'Verbose'                       logical, default false
%
% File format
%   Each charge file may contain one non-comment row per atom in one of
%   the following simple whitespace-delimited formats:
%
%       q
%       local_index   q
%       element       q
%       local_index   element   q
%
%   Comment lines beginning with # or % are ignored.
%
% Notes
%   - molIDs may be a scalar (single cation/anion) or a vector.
%   - In uniform mode, the total charge for each molecule is distributed
%     uniformly over its sites.
%   - In file mode, template-order charges are mapped onto the selected
%     supercell molecule using frame-based atom matching.

    p = inputParser;
    addRequired(p, 'sys', @isstruct);
    addRequired(p, 'molIDs', @(x) isnumeric(x) && ~isempty(x));

    addParameter(p, 'Mode', 'uniform', @(x) ischar(x) || isstring(x));
    addParameter(p, 'TotalCharges', [], @isnumeric);

    addParameter(p, 'ChargeFiles', {}, @(x) iscell(x) || isstring(x) || ischar(x));
    addParameter(p, 'TemplateCoords', {}, @iscell);
    addParameter(p, 'TemplateSpecies', {}, @iscell);

    addParameter(p, 'ReferenceAxis', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 3));
    addParameter(p, 'DistanceTol', 1e-2, @(x) isnumeric(x) && isscalar(x) && x >= 0);

    addParameter(p, 'SetActive', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'DisablePolarizabilityOnCharged', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'ZeroExistingCharges', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));

    parse(p, sys, molIDs, varargin{:});
    opt = p.Results;

    validate_sys(sys);

    molIDs = molIDs(:);
    nMol = numel(molIDs);

    mode = lower(char(string(opt.Mode)));
    if ~ismember(mode, {'uniform','file'})
        error('builder:apply_molecule_charges:BadMode', ...
            'Mode must be ''uniform'' or ''file''.');
    end

    % Default reference axis = supercell b vector
    if isempty(opt.ReferenceAxis)
        refAxis = sys.super_lattice(2, :);
    else
        refAxis = reshape(opt.ReferenceAxis, 1, 3);
    end

    if opt.ZeroExistingCharges
        sys.site_charge(:) = 0;
    end

    if opt.SetActive
        sys = builder.select_active_molecules(sys, molIDs);
    end

    switch mode
        case 'uniform'
            totalCharges = normalize_numeric_per_molecule(opt.TotalCharges, nMol, 'TotalCharges');

            for k = 1:nMol
                molID = molIDs(k);
                idx = builder.site_indices_for_molecule(sys, molID);
                qTot = totalCharges(k);

                if isempty(idx)
                    error('builder:apply_molecule_charges:UnknownMolID', ...
                        'Molecule ID %d was not found in sys.site_mol_id.', molID);
                end

                qSite = qTot / numel(idx);
                sys.site_charge(idx) = qSite;

                if opt.DisablePolarizabilityOnCharged
                    sys.site_is_polarizable(idx) = false;
                end

                if opt.Verbose
                    fprintf('Uniform charge assignment:\n');
                    fprintf('  molecule ID       = %d\n', molID);
                    fprintf('  total charge      = %+g\n', qTot);
                    fprintf('  site count        = %d\n', numel(idx));
                    fprintf('  charge per site   = %+g\n', qSite);
                end
            end

        case 'file'
            chargeFiles = normalize_cellstr_per_molecule(opt.ChargeFiles, nMol, 'ChargeFiles');

            if numel(opt.TemplateCoords) ~= nMol
                error('builder:apply_molecule_charges:BadTemplateCoords', ...
                    'TemplateCoords must have one entry per molecule.');
            end

            if numel(opt.TemplateSpecies) ~= nMol
                error('builder:apply_molecule_charges:BadTemplateSpecies', ...
                    'TemplateSpecies must have one entry per molecule.');
            end

            for k = 1:nMol
                molID = molIDs(k);
                idx = builder.site_indices_for_molecule(sys, molID);

                if isempty(idx)
                    error('builder:apply_molecule_charges:UnknownMolID', ...
                        'Molecule ID %d was not found in sys.site_mol_id.', molID);
                end

                Xtempl = opt.TemplateCoords{k};
                Stempl = cellstr(opt.TemplateSpecies{k}(:));
                Xtarg = sys.site_pos(idx, :);
                Starg = cellstr(sys.site_type(idx));

                qTemplate = read_charge_file(chargeFiles{k});

                if numel(qTemplate) ~= size(Xtempl, 1)
                    error('builder:apply_molecule_charges:ChargeCountMismatch', ...
                        ['Charge file %s contains %d charges, but template molecule %d ', ...
                         'contains %d atoms.'], ...
                        char(chargeFiles{k}), numel(qTemplate), molID, size(Xtempl, 1));
                end

                match = builder.match_molecule_atoms_by_frame( ...
                    Xtempl, Stempl, Xtarg, Starg, ...
                    'ReferenceAxis', refAxis, ...
                    'DistanceTol', opt.DistanceTol, ...
                    'Verbose', false);

                if ~match.is_within_tolerance
                    warning('builder:apply_molecule_charges:FrameMatchTolerance', ...
                        ['Frame-based atom match for molecule %d exceeds tolerance.\n' ...
                         '  max pair distance = %.6e\n' ...
                         '  rms pair distance = %.6e\n' ...
                         '  tolerance         = %.6e'], ...
                        molID, match.max_pair_distance, match.rms_pair_distance, opt.DistanceTol);
                end

                qTarget = zeros(numel(idx), 1);
                for i = 1:numel(qTemplate)
                    j = match.template_to_target(i);
                    qTarget(j) = qTemplate(i);
                end

                sys.site_charge(idx) = qTarget;

                if opt.DisablePolarizabilityOnCharged
                    sys.site_is_polarizable(idx) = false;
                end

                if opt.Verbose
                    fprintf('File-based charge assignment:\n');
                    fprintf('  molecule ID          = %d\n', molID);
                    fprintf('  charge file          = %s\n', char(chargeFiles{k}));
                    fprintf('  atom count           = %d\n', numel(idx));
                    fprintf('  total assigned q     = %+0.10f\n', sum(qTarget));
                    fprintf('  frame-match max dist = %.6e\n', match.max_pair_distance);
                    fprintf('  frame-match rms dist = %.6e\n', match.rms_pair_distance);
                end
            end
    end

    % Small bookkeeping record
    sys.charged_molecule_ids = molIDs(:).';
    sys.charge_assignment_mode = mode;
end


function validate_sys(sys)
    required = {'site_charge', 'site_is_polarizable', 'site_pos', 'site_type', 'site_mol_id', 'super_lattice'};
    for k = 1:numel(required)
        name = required{k};
        if ~isfield(sys, name) || isempty(sys.(name))
            error('builder:apply_molecule_charges:MissingField', ...
                'sys.%s is required and missing/empty.', name);
        end
    end
end


function vals = normalize_numeric_per_molecule(x, nMol, argName)
    if isempty(x)
        error('builder:apply_molecule_charges:MissingNumericInput', ...
            '%s is required.', argName);
    end

    x = x(:);

    if numel(x) ~= nMol
        error('builder:apply_molecule_charges:BadNumericInput', ...
            '%s must have one value per molecule.', argName);
    end

    vals = x;
end


function out = normalize_cellstr_per_molecule(x, nMol, argName)
    if isempty(x)
        error('builder:apply_molecule_charges:MissingCellInput', ...
            '%s is required.', argName);
    end

    if ischar(x) || isstring(x)
        x = cellstr(x);
    end

    if ~iscell(x)
        error('builder:apply_molecule_charges:BadCellInput', ...
            '%s must be a cell array of filenames.', argName);
    end

    if numel(x) ~= nMol
        error('builder:apply_molecule_charges:BadCellInput', ...
            '%s must have one entry per molecule.', argName);
    end

    out = cell(nMol, 1);
    for i = 1:nMol
        out{i} = char(string(x{i}));
    end
end


function q = read_charge_file(filename)
% Read one charge per non-comment line.
% Accepted row formats:
%   q
%   idx q
%   elem q
%   idx elem q

    if ~isfile(filename)
        error('builder:apply_molecule_charges:MissingChargeFile', ...
            'Charge file not found: %s', filename);
    end

    txt = fileread(filename);
    lines = splitlines(string(txt));

    q = zeros(0, 1);

    for i = 1:numel(lines)
        line = strtrim(lines(i));

        if strlength(line) == 0
            continue;
        end

        if startsWith(line, "#") || startsWith(line, "%")
            continue;
        end

        toks = regexp(char(line), '\s+', 'split');
        toks = toks(~cellfun(@isempty, toks));

        switch numel(toks)
            case 1
                qval = str2double(toks{1});

            case 2
                % Either "idx q" or "elem q"
                qval = str2double(toks{2});

            otherwise
                % Use last column as charge
                qval = str2double(toks{end});
        end

        if ~isfinite(qval)
            error('builder:apply_molecule_charges:BadChargeFile', ...
                'Could not parse charge on line %d of %s.', i, filename);
        end

        q(end+1,1) = qval; %#ok<AGROW>
    end
end