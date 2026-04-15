function run_all_nonperiodic_refactor_tests
clc;

fprintf('============================================================\n');
fprintf(' Running nonperiodic refactor regression suite\n');
fprintf('============================================================\n\n');

tests = { ...
    @run_nonperiodic_refactor_test, ...
    @run_nonperiodic_refactor_test_large, ...
    @run_nonperiodic_refactor_test_depolarized, ...
    @run_nonperiodic_refactor_test_initial_mu, ...
    @run_nonperiodic_active_space_energy_test ...
};

names = { ...
    'run_nonperiodic_refactor_test', ...
    'run_nonperiodic_refactor_test_large', ...
    'run_nonperiodic_refactor_test_depolarized', ...
    'run_nonperiodic_refactor_test_initial_mu', ...
    'run_nonperiodic_active_space_energy_test' ...
};

nTests = numel(tests);
times = zeros(nTests, 1);
status = strings(nTests, 1);
messages = strings(nTests, 1);

suiteStart = tic;

for k = 1:nTests
    fprintf('------------------------------------------------------------\n');
    fprintf(' [%d/%d] %s\n', k, nTests, names{k});
    fprintf('------------------------------------------------------------\n');

    tTest = tic;

    try
        tests{k}();
        times(k) = toc(tTest);
        status(k) = "PASS";
        messages(k) = "";
        fprintf('>>> %s: PASS (%.2f s)\n\n', names{k}, times(k));
    catch ME
        times(k) = toc(tTest);
        status(k) = "FAIL";
        messages(k) = string(ME.message);

        fprintf('>>> %s: FAIL (%.2f s)\n', names{k}, times(k));
        fprintf('    %s\n\n', ME.message);

        fprintf('============================================================\n');
        fprintf(' Regression suite FAILED\n');
        fprintf('============================================================\n');
        rethrow(ME);
    end
end

suiteTime = toc(suiteStart);

fprintf('============================================================\n');
fprintf(' Regression suite summary\n');
fprintf('============================================================\n');

for k = 1:nTests
    fprintf(' %-45s  %-4s  %7.2f s\n', names{k}, status(k), times(k));
end

fprintf('\nTotal time: %.2f s\n', suiteTime);
fprintf('ALL TESTS PASSED\n');
end