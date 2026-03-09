
% Calculate gridness score for an autocorrelogram
%
% Calculates a gridness score by expanding a circle around the centre field and
% calculating a correlation value of the expanded circle with it's rotated versions.
% The expansion is done up until the smallest side of the autocorrelogram.
% Can also calculate grid statistics.
%
%  USAGE
%   [score, <stats>] = analyses.gridnessScore(aCorr, <options>)
%   aCorr       A 2D autocorrelogram.
%   <options>   Optional list of property-value pairs (see table below)
%
%   ==============================================================================================
%    Properties         Values
%   ----------------------------------------------------------------------------------------------
%    'threshold'        Normalized threshold value used to search for peaks on the autocorrelogram.
%                       Ranges from 0 to 1, default value is 0.2.
%    'minOrientation'   Value of minimal difference of inner fields orientation (in degrees). If
%                       there are fields that differ in orientation for less than minOrientation,
%                       then only the closest to the centre field are left. Default value is 15.
%   ==============================================================================================
%   score       Gridness score. Ranges from -2 to 2. 2 is more a theoretical bound for a perfect grid.
%               More practical value is around 1.3.
%   stats       If this variable is requested, then it is a structure with the following statistics:
%       spacing         3-element vector with distances from the centre field to neighbour fields.
%       orientation     3-element vector with orientations between the centre field and neighbour fields.
%       ellipse         Ellipse fitted to the grid. Contains the centre, radii and orientation in
%                       radians, stored as [Cx, Cy, Rx, Ry, theta].
%       ellipseTheta    Radius of the ellipse in degrees wrapped in range [0..180].
%
function [gscore, varargout] = gridnessScore(aCorr, varargin)
    nout = max(nargout, 1) - 1;
    inp = inputParser;
    defaultFieldThreshold = 0;
    defaultMinOrientation = 15;
    gridStat.orientation = [];
    gridStat.spacing = [];
    gridStat.ellipse = [];
    gridStat.ellipseTheta = nan;

    % input argument check functions
    checkThreshold = @(x) helpers.isdscalar(x, '>0', '<1');
    checkDScalar = @(x) helpers.isdscalar(x, '>0');

    % fill input parser object
    addRequired(inp, 'aCorr');
    addParameter(inp, 'threshold', defaultFieldThreshold, checkThreshold);
    addParameter(inp, 'minOrientation', defaultMinOrientation, checkDScalar);

    parse(inp, aCorr, varargin{:});

    % get parsed arguments
    fieldThreshold = inp.Results.threshold;
    minOrientation = inp.Results.minOrientation;

    halfSize = ceil(size(aCorr)/2);
    half_height = halfSize(1);
    half_width = halfSize(2);
    aCorrRad = min(halfSize);
    aCorrSize = size(aCorr);

    if aCorrSize(1) == 1 || aCorrSize(2) == 1
        gscore = nan;
        if nout > 0
            varargout{1} = gridStat;
            varargout{2} = 0;
            varargout{3} = nan(6, 2);
            varargout{4} = 0;
        end
        return;
    end

    % contourc is efficient if aCorr is normalized
    maxValue = max(max(aCorr));
    if maxValue ~= 1
        aCorr = aCorr / maxValue;
    end

    cFieldRadius = findCentreRadius(aCorr, fieldThreshold, half_width, half_height);
    % if cFieldRadius == -1
    %     % let's try with increased threshold, this might help
    %     cFieldRadius = findCentreRadius(aCorr, fieldThreshold + 0.1, half_width, half_height);
    % end
    if cFieldRadius == 0 || cFieldRadius == 1 || cFieldRadius == -1
        gscore = nan;
        if nout > 0
            varargout{1} = gridStat;
            varargout{2} = 0;
            varargout{3} = nan(6, 2);
            varargout{4} = 0;
        end
        return;
    end

    % Meshgrid for expanding circle
    [rr, cc] = meshgrid(1:size(aCorr, 2), 1:size(aCorr, 1));

    % Define iteration radius step size for the gridness score
    radSteps = cFieldRadius:aCorrRad;
    radSteps(1) = [];
    numSteps = length(radSteps);

    GNS = zeros(numSteps, 2);
    rotCorr = zeros(1, 5);
    rotAngles_deg = 30*(1:5);

    rotatedCorr = cell(1, length(rotAngles_deg));
    for i = 1:length(rotAngles_deg)
        rotatedCorr{i} = imrotate(aCorr, rotAngles_deg(i), 'bilinear', 'crop');
    end

    mainCircle = sqrt((cc - half_height).^2 + (rr - half_width).^2);
    innerCircle = mainCircle > cFieldRadius;

    % Define expanding ring of autocorrellogram and do x30 correlations
    for i = 1:numSteps
        ind = (innerCircle & (mainCircle < radSteps(i)));
        tempCorr = reshape(aCorr(ind), 1, [])';
        for j = 1:5
            rotatedCircle = reshape(rotatedCorr{j}(ind), 1, [])';
            rotCorr(j) =  corr(tempCorr, rotatedCircle);
        end
        GNS(i, 1) = min(rotCorr([2, 4])) - max(rotCorr([1, 3, 5]));
        GNS(i, 2) = radSteps(i);
    end

    % Find the biggest gridness score and radius
    [gscore, gscoreLoc] = max(GNS(:, 1));
    numGridnessRadii = 3;
    numStep = numSteps - numGridnessRadii;
    if numStep < 1
        numStep = 1;
    end

    if numStep == 1
        gscore = nanmean(GNS(:, 1));
        % radius = nanmean(gridnessArray(:,2));
    else
        meanGridnessArray = zeros(numStep, 1);
        for ii = 1:numStep
            meanGridnessArray(ii) = nanmean(GNS(ii:ii + numGridnessRadii-1, 1));
        end

        [gscore, gInd] = max(meanGridnessArray);
        %gscoreLoc = gridnessArray(gInd+(numGridnessRadii-1)/2, 2);
        gscoreLoc = gInd + (numGridnessRadii-1)/2;
    end

%     [peakLoc, peakMag] = peakfinder(GNS(:, 1));
%     mid = find(peakLoc > 8, 1);
%     if ~isempty(mid)
%         gscore = peakMag(mid);
%         gscoreLoc = mid;
%     else
%         gscore = peakMag(1);
%         gscoreLoc = 1;
%     end

    varargout{4} = radSteps(gscoreLoc);

    % Return if we do not need to calculate grid statistics
    if nout < 1
        return;
    end

    %% Calculate gridness score statistics
    bestCorr = (mainCircle < radSteps(gscoreLoc) * 1.25) .* aCorr;
    regionalMaxMap = imregionalmax(bestCorr, 4);
    se = strel('square', 3);
    im2 = imdilate(regionalMaxMap, se); % dilate map to eliminate fragmentation
    cc = bwconncomp(im2, 8);
    stats = regionprops(cc, 'Centroid');

    if length(stats) < 5
        warning('BNT:numFields', 'Not enough inner fields has been found. Can''t calculate grid properties');

        varargout{1} = gridStat;
        varargout{2} = cFieldRadius;
        varargout{3} = nan(6, 2);
        varargout{4} = radSteps(gscoreLoc);
        return;
    end

    allCoords = [stats(:).Centroid];
    centresOfMass(:, 1) = allCoords(1:2:end);
    centresOfMass(:, 2) = allCoords(2:2:end);

    % Calculate orientation for each field relative to the centre field
    orientation = (atan2(centresOfMass(:, 2) - half_height, centresOfMass(:, 1) - half_width)); % atan2(Y, X)
    peaksToCentre = sqDistance(centresOfMass', [half_width half_height]');
    zeroInd = find(orientation == 0, 1);
    orientation(zeroInd) = []; % remove zero value, so that we do not have a side effect with minOrientation
    stats(zeroInd) = [];
    peaksToCentre(zeroInd) = [];
    centresOfMass(zeroInd, :) = [];

    % filter fields that have similar orientation
    orientDistSq = circ_dist2(orientation);
    closeFields = abs(orientDistSq) < deg2rad(minOrientation);
    [rows, cols] = size(closeFields);
    closeFields(1:(rows+1):rows*cols) = 0; % assign zero to diagonal elements
    closeFields(tril(true(rows))) = 0; % assign zero to lower triangular of a matrix. Matrix is
                                       % symmetric and we do not need these values.
    [rows, cols] = find(closeFields); % find non-empty elements, they correspond to indices of close fields
    if ~isempty(rows)
        indToDelete = zeros(1, length(rows));
        for i = 1:length(rows)
            % fieldPeaks = [fields([rows(i) cols(i)]).peakX; fields([rows(i) cols(i)]).peakY];
            % peaksToCentre = sqDistance(fieldPeaks, [half_width; half_height]);
            if peaksToCentre(rows(i)) > peaksToCentre(cols(i))
                indToDelete(i) = rows(i);
            else
                indToDelete(i) = cols(i);
            end
        end
        indToDelete = unique(indToDelete);
        stats(indToDelete) = [];
        peaksToCentre(indToDelete) = [];

        if length(stats) < 4
            warning('BNT:numFields', 'Not enough inner fields has been found. Can''t calculate grid properties');

            varargout{1} = gridStat;
            varargout{2} = cFieldRadius;
            varargout{3} = nan(6, 2);
            varargout{4} = radSteps(gscoreLoc);
            return;
        end

        allCoords = [stats(:).Centroid];
        clear centresOfMass;
        centresOfMass(:, 1) = allCoords(1:2:end);
        centresOfMass(:, 2) = allCoords(2:2:end);
    end

    % % get fields peak coordinates
    % fieldPeaks = zeros(length(stats), 2);
    % for i = 1:length(stats)
    %     [~, maxInd] = max(bestCorr(stats(i).PixelIdxList));
    %     fieldPeaks(i, :) = stats(i).PixelList(maxInd, :);
    % end
    % % fieldPeaks = [fields(:).peakX; fields(:).peakY]'; %
%     peaksToCentre = sqDistance(centresOfMass', [half_width half_height]');
    [~, sortInd] = sort(peaksToCentre);
    stats = stats(sortInd);
    centresOfMass = centresOfMass(sortInd, :);

    % leave only 6 closest neighbours (if available)
    if length(stats) > 5
%         stats = stats(1:6);
        centresOfMass = centresOfMass(1:6, :);
    else
%         stats = stats(1:end);
        centresOfMass = centresOfMass(1:end, :);
    end

    % centresOfMass = [fields(:).x; fields(:).y]';
    % Calculate orientation for each field relative to the centre field
    orientation = rad2deg(atan2(centresOfMass(:, 2) - half_height, centresOfMass(:, 1) - half_width)); % atan2(Y, X)

    % Calculate distances between centre of masses for each field and the centre field
    spacing = sqrt((centresOfMass(:, 1) - half_width).^2 + (centresOfMass(:, 2) - half_height).^2);

%     % Plot grid polygon points
%     figure, plot.colorMap(bestCorr), hold on;
%     plot(centresOfMass(:, 1), centresOfMass(:, 2), '+k', 'markersize', 8);

    ell = general.fitEllipse(centresOfMass(:, 1), centresOfMass(:, 2));
    ellipseTheta = rad2deg(general.wrap(ell(end)) + pi);
%     drawEllipse(ell, 'linewidth', 2, 'color', [1 1 1]);

    % Determine axes orientation, spacing and deviation
    [~, bBC] = sort(abs(orientation));
    [~, bBC2] = sort(abs(orientation - orientation(bBC(1))));

    % leave only three values, because autocorrelogram is symmetric
    orientation = orientation(bBC2(1:3));
    spacing = spacing(bBC2(1:3));
    [orientation, orientSortInd] = sort(orientation);

    spacing = spacing(orientSortInd);

    gridStat.orientation = orientation;
    gridStat.spacing = spacing;
    gridStat.ellipse = ell;
    gridStat.ellipseTheta = ellipseTheta;

    varargout{1} = gridStat;
    varargout{2} = cFieldRadius;
    varargout{3} = centresOfMass;
    varargout{4} = radSteps(gscoreLoc);
end

function D = sqDistance(X, Y)
    D = bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y);
end

function cFieldRadius = findCentreRadius(aCorr, fieldThreshold, half_width, half_height)
%     figure, contour(aCorr, [fieldThreshold fieldThreshold]), hold on;

    % create a contour plot of fields
    cField = contourc(aCorr, [fieldThreshold fieldThreshold]);
    [~, fLoc] = find(cField(1, :) == fieldThreshold);
    if length(fLoc) > 1;
        fLocMeans = zeros(length(fLoc), 2);
        allFields = cell(length(fLoc), 1);
        for i = 1:length(fLoc)-1 % the last contour is processed differently.
            allFields{i} = cField(:, fLoc(i)+1:fLoc(i+1)-1);
            fLocMeans(i, :) = [mean(allFields{i}(1, :)) mean(allFields{i}(2, :))];
%             text(fLocMeans(i, 1), fLocMeans(i, 2), num2str(i));
        end
        i = i + 1;
        allFields{i} = cField(:, fLoc(i)+1:end);
        fLocMeans(i, :) = [mean(allFields{i}(1, :)) mean(allFields{i}(2, :))];

        % The 'min' approach is a bit faster (~ 0.02s for 1000 calculations), however sqDistance
        % appraoch is somewhat 'more correct'.
        % !!! Min approach gives an error on this data:
        % points = helpers.hexGrid([0 0 50 50], 15);
        % rmap = helpers.gauss2d(points, 10*ones(size(points, 1), 1), [50 50]);
        % aCorr = xcorr2(rmap - mean(rmap(:)));

%         [~, fLocMin] = min(abs(mean([fLocMeans(:, 1) - half_width fLocMeans(:, 2) - half_height], 2)));
%         [~, fLocMin] = min(sqDistance(fLocMeans', [half_width half_height]')); % point should be in format [x y]

        % get all distances and check two minimums of them
        allDistances = sqDistance(fLocMeans', [half_width half_height]'); % point should be in format [x y]
        [~, sortIndices] = sort(allDistances);
        twoMinIndices = sortIndices(1:2);

        if abs(allDistances(twoMinIndices(1)) - allDistances(twoMinIndices(2))) < 1
            % two fields with close middle points. Let's select one with minimum square
            areas = zeros(length(twoMinIndices), 1);
            for i = 1:length(twoMinIndices)
                testedField = allFields{twoMinIndices(i)};
                areas(i) = polyarea(testedField(1, :), testedField(2, :));

                % check that this polygon actually contains the middle point
                % if ~inpolygon(half_width, half_height, testedField(1, :), testedField(2, :))
                %     areas(i) = Inf;
                % end
            end
            % if all(isinf(areas))
            %     cFieldRadius = -1;
            %     return;
            % end
            [~, fLocMin] = min(areas);
            fLocMin = twoMinIndices(fLocMin);
        else
            fLocMin = twoMinIndices(1); % get the first minimum
        end

        centerField = allFields{fLocMin};
    else
        centerField = cField(:, fLoc+1:end);
    end
    cFieldRadius = floor(sqrt(polyarea(centerField(1, :), centerField(2, :))/pi));
end