classdef GridnessScoreTester < matlab.unittest.TestCase
    methods(Test)
        function noArguments(testCase)
            testCase.verifyError(@()analyses.gridnessScore(), 'MATLAB:minrhs');
        end

        function nanMap(testCase)
            map = NaN(50);
            gscore = analyses.gridnessScore(map);

            testCase.verifyEqual(gscore, nan);
        end

        function zeroMap(testCase)
            map = zeros(50);
            gscore = analyses.gridnessScore(map);

            testCase.verifyEqual(gscore, nan);
        end

        function image(testCase)
            % Image from a Matlab bundle
            I = double(rgb2gray(imread('fabric.png')));
            gscore = analyses.gridnessScore(I);

            testCase.verifyEqual(gscore, nan);
        end

        function perfectScore(testCase)
            points = helpers.hexGrid([0 0 50 50], 15);
            rmap = helpers.gauss2d(points, 5*ones(size(points, 1), 1), [50 50]);
            aCorr = analyses.autocorrelation(rmap);

            gscore = analyses.gridnessScore(aCorr);
            gscore = round(gscore * 1000) / 1000;

            testCase.verifyEqual(gscore, 1.399);
        end

        function threeFieldsMap(testCase)
            sigma = [8 4 3];
            mu = [14 11; 5 18; 12 21];
            rmap = helpers.gauss2d(mu, sigma', [25 25]);
            aCorr = analyses.autocorrelation(rmap);

            gscore = analyses.gridnessScore(aCorr);
            gscore = round(gscore * 10000) / 10000;

            testCase.verifyEqual(gscore, -0.1416);
        end

        function rectangularMap(testCase)
            sigma = [8 2]';
            mu = [11 25; 4 5];
            rmap = helpers.gauss2d(mu, sigma, [50 25]);
            aCorr = analyses.autocorrelation(rmap);

            gscore = analyses.gridnessScore(aCorr);
            gscore1 = round(gscore * 10000) / 10000;
            testCase.verifyEqual(gscore1, -0.0022);

            gscore = analyses.gridnessScore(aCorr');
            gscore2 = round(gscore * 10000) / 10000;
            testCase.verifyEqual(gscore2, -0.0022);
            testCase.verifyLessThan(abs(gscore1 - gscore2), 0.05);
        end
    end
end