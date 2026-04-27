function result = select_same_stack_neighbor(desc, shell, varargin)
%SELECT_SAME_STACK_NEIGHBOR Select same-stack complete molecules from a
% reference-relative descriptor table.
%
% result = builder.select_same_stack_neighbor(desc, shell)
% result = builder.select_same_stack_neighbor(desc, shell, ...)
%
% Inputs
%   desc   output of builder.complete_molecule_descriptors_relative_to_reference
%   shell  positive integer shell index
%
% Optional name-value inputs
%   'PerpTol'      scalar, default 1e-3
%   'Direction'    'either' | '+' | '-', default 'either'
%   'ShellTol'     scalar, default 1e-3
%   'Verbose'      logical, default false
%
% Output
%   result struct with fields:
%       .selector
%       .shell
%       .direction
%       .perp_tol
%       .match_table
%       .candidate_table
%       .shell_values

    p = inputParser;
    addRequired(p, 'desc', @isstruct);
    addRequired(p, 'shell', @(x) isnumeric(x) && isscalar(x) && x >= 1 && mod(x,1) == 0);
    addParameter(p, 'PerpTol', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'Direction', 'either', @(x) ischar(x) || isstring(x));
    addParameter(p, 'ShellTol', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
    parse(p, desc, shell, varargin{:});
    opt = p.Results;

    T = validate_desc(desc);

    direction = char(string(opt.Direction));
    if ~ismember(direction, {'either','+','-'})
        error('builder:select_same_stack_neighbor:BadDirection', ...
            'Direction must be ''either'', ''+'', or ''-''.');
    end

    % Same-stack = negligible perpendicular displacement
    candidateMask = T.d_perp <= opt.PerpTol;
    Tcand = T(candidateMask, :);

    if isempty(Tcand)
        error('builder:select_same_stack_neighbor:NoSameStackCandidates', ...
            'No same-stack complete molecules found within PerpTol = %g.', opt.PerpTol);
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
        error('builder:select_same_stack_neighbor:NoDirectionalCandidates', ...
            'No same-stack molecules found for direction ''%s''.', direction);
    end

    absPar = abs(Tcand.d_par);
    [absParSorted, order] = sort(absPar, 'ascend');
    Tcand = Tcand(order, :);

    shellIndex = assign_shells(absParSorted, opt.ShellTol);
    maxShell = max(shellIndex);

    if shell > maxShell
        error('builder:select_same_stack_neighbor:ShellUnavailable', ...
            'Requested shell %d, but only %d same-stack shell(s) are available.', ...
            shell, maxShell);
    end

    matchMask = (shellIndex == shell);
    Tmatch = Tcand(matchMask, :);

    shellValues = zeros(maxShell, 1);
    for s = 1:maxShell
        shellValues(s) = min(absParSorted(shellIndex == s));
    end

    result = struct();
    result.selector = 'same_stack';
    result.shell = shell;
    result.direction = direction;
    result.perp_tol = opt.PerpTol;
    result.match_table = Tmatch;
    result.candidate_table = Tcand;
    result.shell_values = shellValues;

    if opt.Verbose
        fprintf('Same-stack selector summary:\n');
        fprintf('  shell requested       = %d\n', shell);
        fprintf('  direction             = %s\n', direction);
        fprintf('  perp tolerance        = %g\n', opt.PerpTol);
        fprintf('  candidates in stack   = %d\n', height(Tcand));
        fprintf('  available shells      = %d\n', maxShell);
        fprintf('  matched molecules     = %s\n', mat2str(Tmatch.molecule_id(:).'));
    end
end


function T = validate_desc(desc)
    if ~isfield(desc, 'table') || isempty(desc.table) || ~istable(desc.table)
        error('builder:select_same_stack_neighbor:BadDesc', ...
            'desc.table must be a nonempty table.');
    end

    T = desc.table;
    needed = {'molecule_id','d_par','d_perp'};
    for k = 1:numel(needed)
        if ~ismember(needed{k}, T.Properties.VariableNames)
            error('builder:select_same_stack_neighbor:BadDesc', ...
                'desc.table is missing required column ''%s''.', needed{k});
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