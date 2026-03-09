classdef MeanRateTester < tests.TestCaseWithData
    methods
        function this = MeanRateTester()
            this.copyRecordings = false;
        end
    end

    methods(Test)
        function noArguments(testCase)
            testCase.verifyError(@()analyses.meanRate(), 'MATLAB:minrhs');
        end

        function emptyArguments(testCase)
            rate = analyses.meanRate([], []);
            testCase.verifyEqual(rate, 0);
        end

        function values(testCase)
            spikes = [1 2 3];
            pos = ones(15, 3);
            pos(:, 1) = 0:0.02:(0.02*14);
            duration = 0.02 * 15;
            exp = 3/duration;

            rate = analyses.meanRate(spikes, pos);
            testCase.verifyEqual(rate, exp);
        end

        function nonUniformSamplingTime(testCase)
            posFile = fullfile(testCase.DataFolder, 'neuralynx_pos.mat'); % contains variable pos
            load(posFile);

            rate = analyses.meanRate(ones(1000, 1), pos);
            testCase.verifyEqual(roundn(rate, -4), 1.1406);
        end

        function testAxona(testCase)
            posFile = fullfile(testCase.DataFolder, 'axona_pos.mat'); % contains variable pos
            load(posFile);

            rate = analyses.meanRate(ones(1000, 1), pos);
            testCase.verifyEqual(roundn(rate, -4), 0.8319);
        end
    end
end