function result = select_skew_neighbor(desc, shell, varargin)
%SELECT_SKEW_NEIGHBOR Select skew complete molecules from a descriptor table.
%
% A skew neighbor is defined here as a molecule whose normal-angle relative
% to the reference is near a target value (default 72 deg).
%
% result = builder.select_skew_neighbor(desc, shell)
% result = builder.select_skew_neighbor(desc, shell, ...)
%
% Inputs
%   desc   output of builder.complete_molecule_descriptors_relative_to_reference
%          that includes normal_angle_deg
%   shell  positive integer shell index by total distance
%
% Optional name-value inputs
%   'NormalTargetDeg'   scalar, default 72
%   'NormalTolDeg'      scalar, default 15
%   'Direction'         'either' | '+' | '-', default 'either'
%   'ShellTol'          scalar, default 1e-3
%   'Verbose'           logical, default false
%
% Output
%   result struct with fields:
%       .selector
%       .shell
%       .direction
%       .match_table
%       .candidate_table
%       .shell_values

    p = inputParser;
    addRequired(p, 'desc', @isstruct);
    addRequired(p, 'shell', @(x) isnumeric(x) && isscalar(x) && x >= 1 && mod(x,1) == 0);
    addParameter(p, 'NormalTargetDeg', 72, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'NormalTolDeg', 15, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'Direction', 'either', @(x) ischar(x) || isstring(x));
    addParameter(p, 'ShellTol', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
    parse(p, desc, shell, varargin{:});
    opt = p.Results;

    T = validate_desc(desc);

    direction = char(string(opt.Direction));
    if ~ismember(direction, {'either','+','-'})
        error('builder:select_skew_neighbor:BadDirection', ...
            'Direction must be ''either'', ''+'', or ''-''.');
    end

    skewMask = abs(T.normal_angle_deg - opt.NormalTargetDeg) <= opt.NormalTolDeg;
    Tcand = T(skewMask, :);

    if isempty(Tcand)
        error('builder:select_skew_neighbor:NoSkewCandidates', ...
            'No skew candidates found near %.3f +/- %.3f deg.', ...
            opt.NormalTargetDeg, opt.NormalTolDeg);
    end

    switch direction
        case '+'
            Tcand = Tcand(Tcand.d_par > 0, :);
        case '-'
            Tcand = Tcand(Tcand.d_par < 0, :);
        case 'either'
            % keep both
    end

    if isempty(Tcand)
        error('builder:select_skew_neighbor:NoDirectionalCandidates', ...
            'No skew candidates found for direction ''%s''.', direction);
    end

    [distSorted, order] = sort(Tcand.distance, 'ascend');
    Tcand = Tcand(order, :);

    shellIndex = assign_shells(distSorted, opt.ShellTol);
    maxShell = max(shellIndex);

    if shell > maxShell
        error('builder:select_skew_neighbor:ShellUnavailable', ...
            'Requested skew shell %d, but only %d shell(s) are available.', ...
            shell, maxShell);
    end

    matchMask = (shellIndex == shell);
    Tmatch = Tcand(matchMask, :);

    shellValues = zeros(maxShell, 1);
    for s = 1:maxShell
        shellValues(s) = min(distSorted(shellIndex == s));
    end

    result = struct();
    result.selector = 'skew';
    result.shell = shell;
    result.direction = direction;
    result.match_table = Tmatch;
    result.candidate_table = Tcand;
    result.shell_values = shellValues;

    if opt.Verbose
        fprintf('Skew selector summary:\n');
        fprintf('  shell requested       = %d\n', shell);
        fprintf('  direction             = %s\n', direction);
        fprintf('  normal target         = %g deg\n', opt.NormalTargetDeg);
        fprintf('  normal tolerance      = %g deg\n', opt.NormalTolDeg);
        fprintf('  skew candidates       = %d\n', height(Tcand));
        fprintf('  matched molecules     = %s\n', mat2str(Tmatch.molecule_id(:).'));
    end
end


function T = validate_desc(desc)
    if ~isfield(desc, 'table') || isempty(desc.table) || ~istable(desc.table)
        error('builder:select_skew_neighbor:BadDesc', ...
            'desc.table must be a nonempty table.');
    end

    T = desc.table;
    needed = {'molecule_id','distance','d_par','normal_angle_deg'};
    for k = 1:numel(needed)
        if ~ismember(needed{k}, T.Properties.VariableNames)
            error('builder:select_skew_neighbor:BadDesc', ...
                ['desc.table is missing required column ''%s''.\n' ...
                 'Add molecular normals / normal_angle_deg to the descriptor builder first.'], ...
                needed{k});
        end
    end
end


function shellIndex = assign_shells(values, tol)
    n = numel(values);
    shellIndex = zeros(n, 1);

    if n == 0
        return;
    end

    shellCount = 1;
    shellIndex(1) = shellCount;
    rep = values(1);

    for k = 2:n
        if abs(values(k) - rep) <= tol
            shellIndex(k) = shellCount;
        else
            shellCount = shellCount + 1;
            rep = values(k);
            shellIndex(k) = shellCount;
        end
    end
end