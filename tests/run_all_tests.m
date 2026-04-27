function run_all_tests()
%RUN_ALL_TESTS Run all Polarize tests in this folder.
%
% Discovers files named test_*.m in the same directory as this file.
% Each test file should define a function with the same name as the file,
% e.g.
%
%   tests/test_geom_core.m
%
% should start with:
%
%   function test_geom_core()
%
% A test passes if it completes without throwing an error.

fprintf('\n============================================================\n');
fprintf('Polarize test suite\n');
fprintf('============================================================\n');

thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);

% Make sure src/ and tests/ are on the MATLAB path when run from anywhere.
repoRoot = fileparts(testDir);
srcDir = fullfile(repoRoot, 'src');

if exist(srcDir, 'dir')
    addpath(genpath(srcDir));
end

addpath(testDir);

files = dir(fullfile(testDir, '**', 'test_*.m'));

if isempty(files)
    fprintf('\nNo test_*.m files found in:\n  %s\n', testDir);
    return;
end

% Sort for deterministic order.
[~, order] = sort({files.name});
files = files(order);

nPass = 0;
nFail = 0;

failures = struct( ...
    'name', {}, ...
    'message', {}, ...
    'identifier', {}, ...
    'stack', {} ...
);

for i = 1:numel(files)
    [~, testName] = fileparts(files(i).name);

    fprintf('\n[%d/%d] %s\n', i, numel(files), testName);

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
        fprintf('\n%d) %s\n', i, failures(i).name);

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