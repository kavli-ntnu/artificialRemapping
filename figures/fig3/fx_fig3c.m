function fx_fig3c

global hippoGlobe 
CYL

for iDataSet = 1:4
    if iDataSet == 1
        load masterMat_AB1
    elseif iDataSet == 2
        load masterMat_AB2
    elseif iDataSet == 3
        load masterMat_AB3
    elseif iDataSet == 4
        load masterMat_AB3
    end
    
    raw = dataOutput;
    
    %% functional thresholds
    thresh.FR = 7;
    
    %% store names of sessions
    sessions = unique(extract.cols(raw,labels,'session'),'stable');
    numSesh = numel(sessions);
    
    %% fix quality
    badQ = '3';
    offQ = '4';
    if iDataSet == 1 || iDataSet == 2
        raw = cleanUpQualityNovel(raw,labels,sessions,numSesh,badQ,offQ);
    elseif iDataSet == 3 || iDataSet == 4
        raw = cleanUpQualityTilRom(raw,labels,sessions,numSesh,badQ,offQ);
    end
    
    %% only look at excitatory cells
    toRemove = unique(extract.cols(raw,labels,'cell num','session',sessions{1},'mean rate','>=',thresh.FR));
    raw = extract.rows(raw,labels,'remove','cell num',toRemove);
    
    %% stability
    tempWithin = extract.cols(raw,labels,'cell num','cc 1 vs 1','>=',0.5);
    tempAcross = extract.cols(raw,labels,'cell num','session',sessions{1},'spatial info','>',0.67);
    toKeep = [tempWithin; tempAcross];
    toKeep = unique(toKeep);
    raw = extract.rows(raw,labels,'keep','cell num',toKeep);
     
    %% current data
    if iDataSet == 1 || iDataSet == 2
        currentData = raw;
    elseif iDataSet == 3
        ctrlCNO = extract.cols(raw,labels,':','genotype','ctrl','injection','CNO');
        currentData = ctrlCNO;
    elseif iDataSet == 4
        DPsaline = extract.cols(raw,labels,':','genotype','DP','injection','sal');
        currentData = DPsaline;
    end
    
    %% store data in 3 dimensions as: cell, session comparison, experimental group
    clusterNums = unique(extract.cols(currentData,labels,'cell num'),'stable');
    numClusters = length(clusterNums);
    
    predictionError = nan(numClusters,3);
    errorVals = nan(numClusters,2);
    spatialCorr = nan(numClusters,2);
    distStore = nan(numClusters,1);
    
    for iCluster = 1:numClusters
        errorVals(iCluster,1) = str2double(clusterNums{iCluster});
        spatialCorr(iCluster,1) = str2double(clusterNums{iCluster});
        
        qBL = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{1});

        if ~strcmpi(qBL,badQ) && ~strcmpi(qBL,offQ)
            % find number of fields
            mapBL = extract.cols(currentData,labels,'rate map','cell num',clusterNums{iCluster},'session',sessions{1});
            [~,fields] = analyses.placefield(mapBL,'threshold',0.4,'binWidth',hippoGlobe.binWidth,'minBins',10,'minPeak',0.5);
            
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
                
                % plot BL rate map & primary peak
%                 figure;
%                 subplot(121)
%                 hold on
%                 colorMapBRK(mapBL);
%                 plot(BLprimaryPeakX,BLprimaryPeakY,'+','markeredgecolor','k','markersize',15,'linewidth',2.5)
                
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
%                     plot(x(secondIdx),y(secondIdx),'*','markeredgecolor','m','markersize',15,'linewidth',2.5)
                end
                if length(ratesSorted) > 2
                    thirdIdx = find(ratesSorted(3) == meanRates);
%                     plot(x(thirdIdx),y(thirdIdx),'x','markeredgecolor','b','markersize',15,'linewidth',2.5)
                end
                if length(ratesSorted) > 3
                    fourthIdx = find(ratesSorted(4) == meanRates);
%                     plot(x(fourthIdx),y(fourthIdx),'d','markeredgecolor','y','markersize',15,'linewidth',2.5)
                end
                
                % store BL primary peak
                prevPeakX = BLprimaryPeakX;
                prevPeakY = BLprimaryPeakY;
                
                % CNO
                for iSession = 2
                    qCNO = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{iSession});
                    mapCNO = extract.cols(currentData,labels,'rate map','cell num',clusterNums{iCluster},'session',sessions{iSession});
                    peakRateCNO = extract.cols(currentData,labels,'peak rate','cell num',clusterNums{iCluster},'session',sessions{iSession});
                    
                    % randCell = clusterNums{randi(numClusters,1,1)};
                    % qCNO = extract.cols(currentData,labels,'quality','cell num',randCell,'exptHalf',expHalves{1},'session',sessions{iSession});
                    % mapCNO = extract.cols(currentData,labels,'rate map','cell num',randCell,'exptHalf',expHalves{1},'session',sessions{iSession});
                    % peakRateCNO = extract.cols(currentData,labels,'peak rate','cell num',randCell,'exptHalf',expHalves{1},'session',sessions{iSession});
                    
                    [~,CNOfields] = analyses.placefield(mapCNO,'threshold',0.4,'binWidth',hippoGlobe.binWidth,'minBins',10,'minPeak',0.5);
                    
                    if ~isempty(CNOfields) && peakRateCNO > 0.5 && ~strcmpi(qCNO,badQ)
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
                        
                        % plot CNO map & location of primary peaks
%                         subplot(1,2,iSession)
%                         hold on
%                         colorMapBRK(mapCNO);
%                         plot(CNOprimaryPeakX,CNOprimaryPeakY,'+','markeredgecolor','k','markersize',15,'linewidth',2.5)
                        
                        % calculate distance b/t current & previous peak
                        distToCNO = pdist([prevPeakX,prevPeakY; CNOprimaryPeakX,CNOprimaryPeakY]);
                        distStore(iCluster,1) = distToCNO;
%                         text(5,35,sprintf('Distance = %.2f',distToCNO));
                        
                        sc = analyses.spatialCrossCorrelation(mapBL,mapCNO);
                        % for cells that remap
                        if distToCNO > 10 % && sc < 0.26
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
                            errorVals(iCluster,iSession) = nanmin(predictionError(iCluster,:));
%                             text(5,40,sprintf('Error = %.2f',nanmin(predictionError(iCluster,:))));
                            
                            spatialCorr(iCluster,2) = analyses.spatialCrossCorrelation(mapBL,mapCNO);
%                             text(5,45,sprintf('Corr = %.2f',spatialCorr(iCluster,2)));  

                            % peak firing in predicted field
                            minDist = nanmin(predictionError(iCluster,:));
                            idx = find(predictionError(iCluster,:) == minDist);
                            
                            if ~isempty(idx)
                                ratePredicted(iCluster,2) = ratesSorted(idx(1)+1);
                            end
                        end
                    end
                end
            end
        end
    end
    
    if iDataSet == 1
        predictionsNovFix = predictionError;
        errorValsNovFix = errorVals;
        spatialCorrNovFix = spatialCorr;
        distStoreFix = distStore;
    elseif iDataSet == 2
        predictionsNov03 = predictionError;
        errorValsNov03 = errorVals;
        spatialCorrNov03 = spatialCorr;
        distStore03 = distStore;
    elseif iDataSet == 3
        predictionsCtrlCNO = predictionError;
        errorValsCtrlCNO = errorVals;
        spatialCorrCtrlCNO = spatialCorr;
        distStoreCtrl = distStore;
    elseif iDataSet == 4
        predictionsDPsaline = predictionError;
        errorValsDPsaline = errorVals;
        spatialCorrDPsaline = spatialCorr;
        distStoreDP = distStore;
    end
end                   

% pool groups
predictionsTotal = [predictionsNovFix; predictionsNov03; predictionsCtrlCNO; predictionsDPsaline];
errorValsTotal = [errorValsNovFix; errorValsNov03; errorValsCtrlCNO; errorValsDPsaline];
spatialCorrTotal = [spatialCorrNovFix; spatialCorrNov03; spatialCorrCtrlCNO; spatialCorrDPsaline];
distStoreTotal = [distStoreFix; distStore03; distStoreCtrl; distStoreDP];

%%
sum(~isnan(errorValsTotal))
nov = errorValsTotal(:,2) * hippoGlobe.binWidth;
nanmedian(nov)

binRange = 0:2:60;
counts = histcounts(nov,binRange) / sum(~isnan(nov));
figure;
bar(1:2:59,counts,'barwidth',1,'edgecolor','k','facecolor',[0.5 0.5 0.5])
set(gca,'xtickmode','manual','xtick',0:10:60,'box','off','fontsize',14,'fontweight','bold')
xlabel('Prediction error (cm)')
ylabel('Proportion of neurons')
ylim([0 0.16])
