
global hippoGlobe
CYL
startup
PLOT = true;

load fig3a_predictionMaps

numClusters = size(mapStore,1);
predictionError = nan(numClusters,3);
errorVals = nan(numClusters,2);
spatialCorr = nan(numClusters,2);
ratePredicted = nan(numClusters,2);

for iCluster = 1:numClusters
    errorVals(iCluster,1) = (clusterNums{iCluster});
    spatialCorr(iCluster,1) = (clusterNums{iCluster});
    ratePredicted(iCluster,1) = (clusterNums{iCluster});
    
    % find number of fields
    mapBL = mapStore{iCluster,1};
    [fieldMap,fields] = analyses.placefield(mapBL,'threshold',0.4,'binWidth',hippoGlobe.binWidth,'minBins',10,'minPeak',0.5);
    
    % BL
    if ~isempty(fields)
        % find location of BL primary peak
        peaks = [];
        for iFieldNum = 1:size(fields,2)
            peaks(iFieldNum) = fields(1,iFieldNum).peak;
        end
        primary = find(peaks == max(peaks));
        if length(primary) > 1
            primary = primary(1);
        end
        BLprimaryPeakX = fields(1,primary).peakX;
        BLprimaryPeakY = fields(1,primary).peakY;
        
        % find local max from BL rate map
        map = mapBL;
        map(isnan(map)) = 0;
        regmax = imregionalmax(map);
        [y x] = find(regmax == 1);
        
        meanRates = [];
        for i = 1:length(x)
            meanRates(i) = mapBL(y(i),x(i));
        end
        
        % find location of second & third peaks in BL & plot
        ratesSorted = sort(meanRates,'descend');
        
        if length(ratesSorted) > 1
            secondIdx = find(ratesSorted(2) == meanRates);
            fieldCheck(iCluster,1) = fieldMap(x(secondIdx),y(secondIdx));
        else
            secondIdx = nan;
        end
        
        if length(ratesSorted) > 2
            thirdIdx = find(ratesSorted(3) == meanRates);
            fieldCheck(iCluster,2) = fieldMap(x(thirdIdx),y(thirdIdx));
        else
            thirdIdx = nan;
        end
        
        if length(ratesSorted) > 3
            fourthIdx = find(ratesSorted(4) == meanRates);
            fieldCheck(iCluster,3) = fieldMap(x(fourthIdx),y(fourthIdx));
        else
            fourthIdx = nan;
        end
        
        % store BL primary peak
        prevPeakX = BLprimaryPeakX;
        prevPeakY = BLprimaryPeakY;
        
        % plot BL rate map & primary peak
        if PLOT == true
            figure;
            subplot(121)
            hold on
            colorMapBRK(mapBL);
            plot(BLprimaryPeakX,BLprimaryPeakY,'+','markeredgecolor','k','markersize',15,'linewidth',2.5)
            plot(x(secondIdx),y(secondIdx),'*','markeredgecolor','m','markersize',15,'linewidth',2.5)
            plot(x(thirdIdx),y(thirdIdx),'x','markeredgecolor','g','markersize',15,'linewidth',2.5)
            plot(x(fourthIdx),y(fourthIdx),'d','markeredgecolor','y','markersize',15,'linewidth',2.5)
        end
        
        % CNO
        mapCNO = mapStore{iCluster,2};
        peakRateCNO = nanmax(nanmax(mapCNO));
        [~,CNOfields] = analyses.placefield(mapCNO,'threshold',0.4,'binWidth',hippoGlobe.binWidth,'minBins',10,'minPeak',0.5);
        
        if ~isempty(CNOfields) && peakRateCNO > 0.5
            % find location of peak for CNO session
            Peak = [];
            for iFieldNum = 1:size(CNOfields,2)
                Peak(iFieldNum) = CNOfields(1,iFieldNum).peak;
            end
            primary = find(Peak == max(Peak));
            if length(primary) > 1
                primary = primary(1);
            end
            CNOprimaryPeakX = CNOfields(1,primary).peakX;
            CNOprimaryPeakY = CNOfields(1,primary).peakY;
            
            % calculate distance b/t current & previous peak
            distToCNO = pdist([prevPeakX,prevPeakY; CNOprimaryPeakX,CNOprimaryPeakY]);
            distStore(iCluster,1) = distToCNO;
            sc = analyses.spatialCrossCorrelation(mapBL,mapCNO);
            
            % distance b/t current & previous second peak
            if length(ratesSorted) > 1
                predictionError(iCluster,1) = pdist([x(secondIdx),y(secondIdx); CNOprimaryPeakX,CNOprimaryPeakY]);
            end
            % distance b/t current & previous third peak
            if length(ratesSorted) > 2
                predictionError(iCluster,2) = pdist([x(thirdIdx),y(thirdIdx); CNOprimaryPeakX,CNOprimaryPeakY]);
            end
            % distance b/t current & previous fourth peak
            if length(ratesSorted) > 3
                predictionError(iCluster,3) = pdist([x(fourthIdx),y(fourthIdx); CNOprimaryPeakX,CNOprimaryPeakY]);
            end
            
            % find minimum distance
            errorVals(iCluster,2) = nanmin(predictionError(iCluster,:));
            spatialCorr(iCluster,2) = analyses.spatialCrossCorrelation(mapBL,mapCNO);
            
            % peak firing in predicted field
            minDist = nanmin(predictionError(iCluster,:));
            idx = find(predictionError(iCluster,:) == minDist);
            if ~isempty(idx)
                ratePredicted(iCluster,2) = ratesSorted(idx(1)+1);
            end
            
            % plot CNO map & location of primary peaks
            if PLOT == true
                subplot(1,2,2)
                hold on
                colorMapBRK(mapCNO);
                plot(CNOprimaryPeakX,CNOprimaryPeakY,'+','markeredgecolor','k','markersize',15,'linewidth',2.5)
                text(5,40,sprintf('Error = %.2f',nanmin(predictionError(iCluster,:))*2))
            end
        end
    end
end


