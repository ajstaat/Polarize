function result = select_side_stack_neighbor(desc, shell, member, varargin)
%SELECT_SIDE_STACK_NEIGHBOR Select side-stack complete molecules from a
% reference-relative descriptor table.
%
% Side-stack molecules are off-stack molecules grouped by transverse
% displacement family in the plane perpendicular to the chosen stack axis.
%
% result = builder.select_side_stack_neighbor(desc, shell, member)
% result = builder.select_side_stack_neighbor(desc, shell, member, ...)
%
% Inputs
%   desc    output of builder.complete_molecule_descriptors_relative_to_reference
%   shell   positive integer transverse shell index
%   member  positive integer longitudinal member index within that shell
%
% Optional name-value inputs
%   'SameStackPerpTol'  scalar, default 1e-3
%   'Direction'         'either' | '+' | '-', default 'either'
%   'ShellTol'          scalar, default 1e-3
%   'VectorGroupTol'    scalar, default 1e-3
%   'Verbose'           logical, default false
%
% Output
%   result struct with fields:
%       .selector
%       .shell
%       .member
%       .direction
%       .match_table
%       .candidate_table
%       .family_table

    p = inputParser;
    addRequired(p, 'desc', @isstruct);
    addRequired(p, 'shell', @(x) isnumeric(x) && isscalar(x) && x >= 1 && mod(x,1) == 0);
    addRequired(p, 'member', @(x) isnumeric(x) && isscalar(x) && x >= 1 && mod(x,1) == 0);
    addParameter(p, 'SameStackPerpTol', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'Direction', 'either', @(x) ischar(x) || isstring(x));
    addParameter(p, 'ShellTol', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'VectorGroupTol', 1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
    parse(p, desc, shell, member, varargin{:});
    opt = p.Results;

    T = validate_desc(desc);

    direction = char(string(opt.Direction));
    if ~ismember(direction, {'either','+','-'})
        error('builder:select_side_stack_neighbor:BadDirection', ...
            'Direction must be ''either'', ''+'', or ''-''.');
    end

    % Side-stack = not same-stack
    sideMask = T.d_perp > opt.SameStackPerpTol;
    Tcand = T(sideMask, :);

    if isempty(Tcand)
        error('builder:select_side_stack_neighbor:NoSideStackCandidates', ...
            'No side-stack candidates found.');
    end

    % Group by transverse displacement family
    V = [Tcand.d_perp_x, Tcand.d_perp_y, Tcand.d_perp_z];
    familyID = group_vector_families(V, opt.VectorGroupTol);

    Tcand.family_id = familyID;

    uniqueFamilies = unique(familyID(:)).';
    nFam = numel(uniqueFamilies);
    famPerpNorm = zeros(nFam, 1);

    for f = 1:nFam
        fam = uniqueFamilies(f);
        row = find(Tcand.family_id == fam, 1, 'first');
        famPerpNorm(f) = Tcand.d_perp(row);
    end

    famShell = assign_shells(famPerpNorm, opt.ShellTol);

    % Store family shell on rows
    Tcand.transverse_shell = zeros(height(Tcand), 1);
    for f = 1:nFam
        fam = uniqueFamilies(f);
        Tcand.transverse_shell(Tcand.family_id == fam) = famShell(f);
    end

    maxShell = max(Tcand.transverse_shell);
    if shell > maxShell
        error('builder:select_side_stack_neighbor:ShellUnavailable', ...
            'Requested side-stack shell %d, but only %d shell(s) are available.', ...
            shell, maxShell);
    end

    Tshell = Tcand(Tcand.transverse_shell == shell, :);

    switch direction
        case '+'
            Tshell = Tshell(Tshell.d_par > 0, :);
        case '-'
            Tshell = Tshell(Tshell.d_par < 0, :);
        case 'either'
            % keep both
    end

    if isempty(Tshell)
        error('builder:select_side_stack_neighbor:NoDirectionalCandidates', ...
            'No side-stack candidates found in shell %d for direction ''%s''.', ...
            shell, direction);
    end

    % Member index is assigned within this shell by |d_par|
    absPar = abs(Tshell.d_par);
    [absParSorted, order] = sort(absPar, 'ascend');
    Tshell = Tshell(order, :);

    memberIndex = assign_shells(absParSorted, opt.ShellTol);
    maxMember = max(memberIndex);

    if member > maxMember
        error('builder:select_side_stack_neighbor:MemberUnavailable', ...
            'Requested member %d in side-stack shell %d, but only %d member shell(s) are available.', ...
            member, shell, maxMember);
    end

    matchMask = (memberIndex == member);
    Tmatch = Tshell(matchMask, :);

    familyTable = table(uniqueFamilies(:), famShell(:), famPerpNorm(:), ...
        'VariableNames', {'family_id','transverse_shell','family_d_perp'});

    result = struct();
    result.selector = 'side_stack';
    result.shell = shell;
    result.member = member;
    result.direction = direction;
    result.match_table = Tmatch;
    result.candidate_table = Tcand;
    result.family_table = familyTable;

    if opt.Verbose
        fprintf('Side-stack selector summary:\n');
        fprintf('  shell requested        = %d\n', shell);
        fprintf('  member requested       = %d\n', member);
        fprintf('  direction              = %s\n', direction);
        fprintf('  side-stack candidates  = %d\n', height(Tcand));
        fprintf('  families found         = %d\n', nFam);
        fprintf('  matched molecules      = %s\n', mat2str(Tmatch.molecule_id(:).'));
    end
end


function T = validate_desc(desc)
    if ~isfield(desc, 'table') || isempty(desc.table) || ~istable(desc.table)
        error('builder:select_side_stack_neighbor:BadDesc', ...
            'desc.table must be a nonempty table.');
    end

    T = desc.table;
    needed = {'molecule_id','d_par','d_perp','d_perp_x','d_perp_y','d_perp_z'};
    for k = 1:numel(needed)
        if ~ismember(needed{k}, T.Properties.VariableNames)
            error('builder:select_side_stack_neighbor:BadDesc', ...
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


function familyID = group_vector_families(V, tol)
    n = size(V, 1);
    familyID = zeros(n, 1);

    if n == 0
        return;
    end

    nFam = 0;
    reps = zeros(0, 3);

    for i = 1:n
        assigned = false;

        for f = 1:nFam
            if norm(V(i,:) - reps(f,:)) <= tol
                familyID(i) = f;
                assigned = true;
                break;
            end
        end

        if ~assigned
            nFam = nFam + 1;
            reps(nFam,:) = V(i,:);
            familyID(i) = nFam;
        end
    end
end