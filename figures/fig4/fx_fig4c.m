load fig4b

spatialCorrAll = nan(200*numSims,1);
for i = 1:200*numSims
    if iBL(i) || iCNO(i)
        mapBL = general.smooth(pBL(:,:,i),6.67);
        mapCNO = general.smooth(pCNO(:,:,i),6.67);
        spatialCorrAll(i,1) = analyses.spatialCrossCorrelation(mapBL,mapCNO);
    end
end

nanmedian(spatialCorrAll)
sum(~isnan(spatialCorrAll))

%%
global hippoGlobe
BIGBOX 
RAND = false; % turn on for running shuffle
PLOT = false;
if PLOT
    mkdir('maps')
end

predictionError = nan(5000,3);
errorVals = nan(5000,1);
spatialCorr = nan(5000,1);
distToCNOstore = nan(5000,1);

for i = 1:5000
    if (iBL(i) || iCNO(i))
        sc = nan;
        mapBL = general.smooth(pBL(:,:,i),6.67);
%         mapBL = L1.p3R(:,:,i);
        [~,fields] = analyses.placefield(mapBL,'threshold',0.1,'minPeak',0.5);
        
        if ~isempty(fields)
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
            if RAND == true
                % random control
                flag = true;
                while flag
                    randIdx = randi(5000);
                    if iCNO(randIdx)
                        flag = false;
                    end
                end
                mapCNO = general.smooth(pCNO(:,:,randIdx),6.67);
            else
                mapCNO = general.smooth(pCNO(:,:,i),6.67);
            end
            
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
                distToCNOstore(i,1) = distToCNO;
                % sc = analyses.spatialCrossCorrelation(mapBL,mapCNO);
                
                if PLOT
                    subplot(1,2,2)
                    hold on
                    colorMapBRK(mapCNO,'clrmap','jet');
                    plot(CNOprimaryPeakX,CNOprimaryPeakY,'+','markeredgecolor','k','markersize',15,'linewidth',2.5)
                end
                
                if distToCNO > 20 % sc < 0.26 % 
                    
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

                end
            end
        end
    end
end

%%
figure;
binRange = 0:3.5:100;
centers = binRange(1:end-1) + diff(binRange)/2;
counts = histcounts(errorVals,binRange) / sum(~isnan(errorVals));
bar(centers,counts,'barwidth',1)
xticks(20:20:100)
% set(gca,'xtickmode','manual','xtick',0:10:100,'box','off','fontsize',14,'fontweight','bold')
xlabel('Prediction offset (cm)')
ylabel('Proportion of neurons')
ylim([0 0.15])

sum(~isnan(errorVals))
nanmedian(errorVals)
