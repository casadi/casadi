classdef ExamplesTests < matlab.unittest.TestCase
    % Unit test suite for CasADi MATLAB examples
    % This test class calls the example .m scripts
    %
    % Excluded files:
    % - MyCallback.m (class definition, not a runnable script)
    %
    % Disabled tests:
    % - testSysidGaussNewton (currently disabled, should be fixed and re-enabled)

    methods(Test, TestTags = {'examples'})
        function testAccessingMxAlgorithm(testCase)
            % Run accessing_mx_algorithm example
            accessing_mx_algorithm
        end

        function testAccessingSxAlgorithm(testCase)
            % Run accessing_sx_algorithm example
            accessing_sx_algorithm
        end

        function testCallback(testCase)
            % Run callback example
            callback
        end

        function testDaebuilder(testCase)
            % Run daebuilder example
            daebuilder
        end

        function testDirectCollocation(testCase)
            % Run direct_collocation example
            direct_collocation
        end

        function testDirectCollocationOpti(testCase)
            % Run direct_collocation_opti example
            direct_collocation_opti
        end

        function testDirectMultipleShooting(testCase)
            % Run direct_multiple_shooting example
            direct_multiple_shooting
        end

        function testDirectSingleShooting(testCase)
            % Run direct_single_shooting example
            direct_single_shooting
        end

        function testLotkaVolterraMinlp(testCase)
            % Run lotka_volterra_minlp example
            lotka_volterra_minlp
        end

        function testRaceCar(testCase)
            % Run race_car example
            race_car
        end

        function testRosenbrock(testCase)
            % Run rosenbrock example
            rosenbrock
        end

        function testSensitivityAnalysis(testCase)
            % Run sensitivity_analysis example
            sensitivity_analysis
        end

        function testSysid(testCase)
            % Run sysid example (not on Windows)
            testCase.assumeFalse(ispc, 'Sysid test is skipped on Windows');
            sysid
        end

        function testSysidGaussNewton(testCase)
            % Run sysid_gauss_newton example (currently disabled)
            testCase.assumeFalse(true, 'sysid_gauss_newton test is currently disabled - needs to be fixed');
            sysid_gauss_newton
        end

        function testVdpIndirectMultipleShooting(testCase)
            % Run vdp_indirect_multiple_shooting example
            vdp_indirect_multiple_shooting
        end
    end
end
