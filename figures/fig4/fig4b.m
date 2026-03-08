load fig4b

%%
cnt = 1;
figure;
selectedCells = [449 466 124 262];
for i = 1:size(selectedCells,2)
    cellNum = selectedCells(i);
    
    mapBL = pBL(:,:,cellNum);
    subplot(4,2,cnt)
    colorMapBRK(general.smooth(mapBL,6.6667));
    title(sprintf('%d',cellNum))
    cnt = cnt + 1;
    
    mapCNO = pCNO(:,:,cellNum);
    subplot(4,2,cnt)
    colorMapBRK(general.smooth(mapCNO,6.6667));
    cnt = cnt + 1;
    
    sc(i,1) = analyses.spatialCrossCorrelation(general.smooth(mapBL,6.6667),general.smooth(mapCNO,6.6667));
end

%%
PLOT = true;

for i = 1:size(selectedCells,2)
    cellNum = selectedCells(i);
    
    mapBL = pBL(:,:,cellNum);
    mapBL = general.smooth(mapBL,6.6667);
    [~,fields] = analyses.placefield(mapBL,'threshold',0.1,'minPeak',0.5);

    if ~isempty(fields)
        fieldSize(i) = nanmax([fields.size]);
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
        
        if PLOT
            figure
            subplot(121)
            hold on
            colorMapBRK(mapBL,'clrmap','jet');
            plot(BLprimaryPeakX,BLprimaryPeakY,'+','markeredgecolor','k','markersize',15,'linewidth',2.5)
        end
        
        % find local max from BL rate map
        map = mapBL;
        map(isnan(map)) = 0;
        regmax = imregionalmax(map);
        [y x] = find(regmax == 1);
        
        meanRates = [];
        for j = 1:length(x)
            meanRates(j) = mapBL(y(j),x(j));
        end
        
        % find location of second & third peaks in BL & plot
        ratesSorted = sort(meanRates,'descend');
        if length(ratesSorted) > 1
            secondIdx = find(ratesSorted(2) == meanRates);
            if PLOT
                plot(x(secondIdx),y(secondIdx),'*','markeredgecolor','m','markersize',15,'linewidth',2.5)
            end
        end
        if length(ratesSorted) > 2
            thirdIdx = find(ratesSorted(3) == meanRates);
            if PLOT
                plot(x(thirdIdx),y(thirdIdx),'x','markeredgecolor','b','markersize',15,'linewidth',2.5)
            end
        end
        if length(ratesSorted) > 3
            fourthIdx = find(ratesSorted(4) == meanRates);
            if PLOT
                plot(x(fourthIdx),y(fourthIdx),'d','markeredgecolor','y','markersize',15,'linewidth',2.5)
            end
        end
        
        % store BL primary peak
        prevMap = mapBL;
        prevPeakX = BLprimaryPeakX;
        prevPeakY = BLprimaryPeakY;
        
        % 2nd map
        mapCNO = pCNO(:,:,cellNum);
        mapCNO = general.smooth(mapCNO,6.6667);
        [~,CNOfields] = analyses.placefield(mapCNO,'threshold',0.1,'minPeak',0.5);
        
        if ~isempty(CNOfields)
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
            sc = analyses.spatialCrossCorrelation(mapBL,mapCNO);
            
            if PLOT
                subplot(1,2,2)
                hold on
                colorMapBRK(mapCNO,'clrmap','jet');
                plot(CNOprimaryPeakX,CNOprimaryPeakY,'+','markeredgecolor','k','markersize',15,'linewidth',2.5)
            end
            
            if distToCNO > 20
                if length(ratesSorted) > 1
                    predictionError(i,1) = pdist([x(secondIdx),y(secondIdx); CNOprimaryPeakX,CNOprimaryPeakY]);
                end
                % distance b/t current & previous third peak
                if length(ratesSorted) > 2
                    predictionError(i,2) = pdist([x(thirdIdx),y(thirdIdx); CNOprimaryPeakX,CNOprimaryPeakY]);
                end
                % distance b/t current & previous fourth peak
                if length(ratesSorted) > 3
                    predictionError(i,3) = pdist([x(fourthIdx),y(fourthIdx); CNOprimaryPeakX,CNOprimaryPeakY]);
                end
                
                % find minimum distance
                errorVals(i) = nanmin(predictionError(i,:));
                if PLOT
                    text(50,-15,sprintf('Error = %.2f',nanmin(predictionError(i,:))));
                end
                spatialCorr(i) = analyses.spatialCrossCorrelation(mapBL,mapCNO);   
                minVal = nanmin(predictionError(i,:));      
            end
        end
    end
end
