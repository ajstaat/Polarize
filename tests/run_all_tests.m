function run_all_tests(varargin)
%RUN_ALL_TESTS Run Polarize tests with optional filtering.
%
% Usage
%   run_all_tests
%   run_all_tests('builder')
%   run_all_tests({'geom','io'})
%   run_all_tests('Only', {'builder','io'})
%   run_all_tests('Exclude', 'viz')
%   run_all_tests('Pattern', 'test_builder_*')
%
% Notes
%   - Tests are discovered recursively under this tests/ folder.
%   - Test files should be named test_*.m.
%   - Each test file should define a function with the same name as the file.
%
% Folder filters match path parts under tests/. For example:
%
%   tests/builder/test_builder_make_crystal_system.m
%
% is matched by 'builder'.
%
% If your tests are still loose directly under tests/, folder filters will
% only work after you move them into subfolders.

fprintf('\n============================================================\n');
fprintf('Polarize test suite\n');
fprintf('============================================================\n');

[thisFile, testDir, repoRoot, srcDir] = local_paths();

if exist(srcDir, 'dir')
    addpath(genpath(srcDir));
end

addpath(genpath(testDir));

opt = local_parse_inputs(varargin{:});

files = dir(fullfile(testDir, '**', opt.Pattern));
files = files(~[files.isdir]);

% Exclude this runner if pattern is broad.
files = files(~strcmp({files.name}, [mfilename '.m']));

files = local_apply_filters(files, testDir, opt);

if isempty(files)
    fprintf('\nNo matching tests found.\n');
    fprintf('  testDir = %s\n', testDir);
    fprintf('  Pattern = %s\n', opt.Pattern);

    if ~isempty(opt.Only)
        fprintf('  Only    = %s\n', strjoin(opt.Only, ', '));
    end

    if ~isempty(opt.Exclude)
        fprintf('  Exclude = %s\n', strjoin(opt.Exclude, ', '));
    end

    return;
end

% Sort deterministically by relative path.
relPaths = local_relative_paths(files, testDir);
[relPaths, order] = sort(relPaths);
files = files(order);

fprintf('\nDiscovered %d test(s).\n', numel(files));

if ~isempty(opt.Only)
    fprintf('Only folders: %s\n', strjoin(opt.Only, ', '));
end

if ~isempty(opt.Exclude)
    fprintf('Excluded folders: %s\n', strjoin(opt.Exclude, ', '));
end

fprintf('\n');

nPass = 0;
nFail = 0;

failures = struct( ...
    'name', {}, ...
    'relpath', {}, ...
    'message', {}, ...
    'identifier', {}, ...
    'stack', {} ...
);

for i = 1:numel(files)
    [~, testName] = fileparts(files(i).name);
    relPath = relPaths{i};

    fprintf('[%d/%d] %s\n', i, numel(files), relPath);

    t = tic;

    try
        feval(testName);
        elapsed = toc(t);

        fprintf('  passed in %.3f s\n', elapsed);
        nPass = nPass + 1;

    catch ME
        elapsed = toc(t);

        fprintf('  FAILED in %.3f s\n', elapsed);
        fprintf('  %s\n', ME.message);

        nFail = nFail + 1;

        failures(end+1).name = testName; %#ok<AGROW>
        failures(end).relpath = relPath;
        failures(end).message = ME.message;
        failures(end).identifier = ME.identifier;
        failures(end).stack = ME.stack;
    end
end

fprintf('\n============================================================\n');
fprintf('Test summary: %d passed, %d failed, %d total\n', ...
    nPass, nFail, numel(files));
fprintf('============================================================\n');

if nFail > 0
    fprintf('\nFailures:\n');

    for i = 1:numel(failures)
        fprintf('\n%d) %s\n', i, failures(i).relpath);

        if ~isempty(failures(i).identifier)
            fprintf('   identifier: %s\n', failures(i).identifier);
        end

        fprintf('   message: %s\n', failures(i).message);

        if ~isempty(failures(i).stack)
            top = failures(i).stack(1);
            fprintf('   location: %s line %d\n', top.name, top.line);
        end
    end

    error('Polarize:TestsFailed', ...
        '%d test(s) failed. See output above.', nFail);
end

end

% =========================================================================
% Local helpers
% =========================================================================

function [thisFile, testDir, repoRoot, srcDir] = local_paths()
thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);
repoRoot = fileparts(testDir);
srcDir = fullfile(repoRoot, 'src');
end

function opt = local_parse_inputs(varargin)

opt = struct();
opt.Only = {};
opt.Exclude = {};
opt.Pattern = 'test_*.m';

if nargin == 0
    return;
end

% Convenience:
%   run_all_tests('builder')
%   run_all_tests({'geom','io'})
if nargin == 1
    arg = varargin{1};

    if ischar(arg) || isstring(arg) || iscellstr(arg)
        opt.Only = local_to_cellstr(arg);
        return;
    end
end

k = 1;

while k <= nargin
    name = varargin{k};

    if ~(ischar(name) || isstring(name))
        error('run_all_tests:BadInput', ...
            'Expected name-value inputs or a folder filter.');
    end

    name = lower(char(string(name)));

    if k == nargin
        error('run_all_tests:MissingValue', ...
            'Missing value for option "%s".', name);
    end

    value = varargin{k+1};

    switch name
        case {'only', 'include', 'folder', 'folders'}
            opt.Only = local_to_cellstr(value);

        case {'exclude', 'skip'}
            opt.Exclude = local_to_cellstr(value);

        case {'pattern', 'filepattern'}
            if ~(ischar(value) || isstring(value))
                error('run_all_tests:BadPattern', ...
                    'Pattern must be a character vector or string scalar.');
            end
            opt.Pattern = char(string(value));

        otherwise
            error('run_all_tests:UnknownOption', ...
                'Unknown option "%s".', name);
    end

    k = k + 2;
end

end

function files = local_apply_filters(files, testDir, opt)

if isempty(files)
    return;
end

relPaths = local_relative_paths(files, testDir);

if ~isempty(opt.Only)
    keep = false(numel(files), 1);

    for i = 1:numel(files)
        parts = local_path_parts(relPaths{i});

        for j = 1:numel(opt.Only)
            if any(strcmpi(parts, opt.Only{j}))
                keep(i) = true;
                break;
            end
        end
    end

    files = files(keep);
    relPaths = relPaths(keep);
end

if ~isempty(opt.Exclude)
    keep = true(numel(files), 1);

    for i = 1:numel(files)
        parts = local_path_parts(relPaths{i});

        for j = 1:numel(opt.Exclude)
            if any(strcmpi(parts, opt.Exclude{j}))
                keep(i) = false;
                break;
            end
        end
    end

    files = files(keep);
end

end

function relPaths = local_relative_paths(files, testDir)

relPaths = cell(numel(files), 1);

for i = 1:numel(files)
    fullPath = fullfile(files(i).folder, files(i).name);

    if startsWith(fullPath, [testDir filesep])
        relPath = extractAfter(fullPath, strlength(testDir) + 1);
    else
        relPath = files(i).name;
    end

    relPaths{i} = char(strrep(relPath, filesep, '/'));
end

end

function parts = local_path_parts(relPath)

relPath = char(string(relPath));
relPath = strrep(relPath, '\', '/');

parts = split(string(relPath), '/');
parts = cellstr(parts(:));

% Drop filename.
if ~isempty(parts)
    parts = parts(1:end-1);
end

end

function c = local_to_cellstr(x)

if ischar(x)
    c = {x};

elseif isstring(x)
    c = cellstr(x(:));

elseif iscellstr(x)
    c = x(:);

else
    error('run_all_tests:BadFilter', ...
        'Filter must be a string, character vector, or cell array of character vectors.');
end

c = cellfun(@(s) char(lower(strtrim(string(s)))), c, 'UniformOutput', false);
c = c(~cellfun(@isempty, c));

end