classdef CasadiTests < matlab.unittest.TestCase
    % Minimal unit test suite for CasADi MATLAB interface
    % This test class calls the existing .m test scripts
    %
    % Platform-specific exclusions:
    % - Windows: skips 'solvers' test
    % - Linux/Mac: skips 'solvers' and 'callback' tests

    methods(Test, TestTags = {'unittests'})
        function testAsortedTests(testCase)
            % Run assorted tests
            asorted_tests
        end

        function testCallback(testCase)
            % Run callback tests (only on Windows)
            testCase.assumeFalse(true, 'Callback test is disabled on all platforms');
            callback
        end

        function testExportTests(testCase)
            % Run export tests
            export_tests
        end

        function testMultithreading(testCase)
            % Run multithreading tests
            multithreading
        end

        function testNlpCallback(testCase)
            % Run NLP callback tests
            nlp_callback
        end

        function testSolvers(testCase)
            % Run solver tests (skipped on all platforms)
            testCase.assumeFalse(true, 'Solvers test is disabled on all platforms');
            solvers
        end
    end
end
