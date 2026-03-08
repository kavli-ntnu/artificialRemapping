function fx_fig3b_3c_3d

global hippoGlobe
CYL
RAND = false; % turn on to generate shuffle
PLOT = false;

% load data
for iDataSet = 1:4
    if iDataSet == 1
        load masterMat_HP1
    elseif iDataSet == 2
        load masterMat_HP2
    elseif iDataSet == 3
        load masterMat_HP3
    elseif iDataSet == 4
        load masterMat_HP4
    end
    raw = dataOutput;
    
    %% store names of groups and sessions
    sessions = unique(extract.cols(raw,labels,'session'),'stable');
    numSesh = numel(sessions);
    
    %% functional thresholds
    thresh.FR = 7;
    thresh.spatInfo = 0.67;
    
    %% fix quality
    badQ = '3';
    offQ = '4';
    if iDataSet == 1
        % switch numbers to strings so you can use functions like 'unique'
        for iRow = 1:size(raw,1)
            raw{iRow,strcmpi('cell num',labels)} = num2str(raw{iRow,strcmpi('cell num',labels)});
            raw{iRow,strcmpi('quality',labels)} = num2str(raw{iRow,strcmpi('quality',labels)});
        end
        
        raw = cleanUpQualityExtended(raw,labels,sessions,numSesh,badQ,offQ);
    elseif iDataSet == 2
        for iRow = 1:size(raw,1)
            raw{iRow,strcmpi('cell num',labels)} = num2str(raw{iRow,strcmpi('cell num',labels)});
            raw{iRow,strcmpi('quality',labels)} = num2str(raw{iRow,strcmpi('quality',labels)});
            raw{iRow,strcmpi('exptHalf',labels)} = num2str(raw{iRow,strcmpi('exptHalf',labels)});
        end
        
        injectionNames = unique(extract.cols(raw,labels,'injection'),'sorted');
        expHalves = unique(extract.cols(raw,labels,'exptHalf'),'stable');
        numExpHalves = numel(expHalves);
        
        raw = cleanUpQualityCML(raw,labels,expHalves,numExpHalves,sessions,numSesh,badQ,offQ);
        raw = extract.cols(raw,labels,':','exptHalf',expHalves{1});
    elseif iDataSet == 3
        raw = cleanUpQualityExtended(raw,labels,sessions,numSesh,badQ,offQ);
    elseif iDataSet == 4
        raw = cleanUpQualityExtended(raw,labels,sessions,numSesh,badQ,offQ);
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
    
    %% subset by group in case you want to plot something now
    if iDataSet == 1
        exp1x = extract.cols(raw,labels,':','group','exp','dose','1x');
    elseif iDataSet == 2
        exp1x = extract.cols(raw,labels,':','group','exp','injection','1x');
    elseif iDataSet == 3
        exp1x = extract.cols(raw,labels,':','genotype','exp','injection','CNO');
    elseif iDataSet == 4
        exp1x = extract.cols(raw,labels,':','genotype','exp','injection','CNO');
    end
    
    %% store data in 3 dimensions as: cell, session comparison, experimental group
    currentData = exp1x;
    clusterNums = unique(extract.cols(currentData,labels,'cell num'),'stable');
    numClusters = length(clusterNums);
    
    predictionError = nan(numClusters,3);
    errorVals = nan(numClusters,6);
    spatialCorr = nan(numClusters,6);
    ratePredicted = nan(numClusters,6);
    distStore = nan(numClusters,1);
    
    for iCluster = 1:numClusters
        
        errorVals(iCluster,1) = str2double(clusterNums{iCluster});
        spatialCorr(iCluster,1) = str2double(clusterNums{iCluster});
        ratePredicted(iCluster,1) = str2double(clusterNums{iCluster});
        
        qBL = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{1});
        calc = nan;
        
        if ~strcmpi(qBL,badQ) && ~strcmpi(qBL,offQ)
            
            % find number of fields
            mapBL = extract.cols(currentData,labels,'rate map','cell num',clusterNums{iCluster},'session',sessions{1});
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
                else
                    secondIdx = nan;
                end
                
                if length(ratesSorted) > 2
                    thirdIdx = find(ratesSorted(3) == meanRates);
                else
                    thirdIdx = nan;
                end
                
                if length(ratesSorted) > 3
                    fourthIdx = find(ratesSorted(4) == meanRates);
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
                if iDataSet == 1 || iDataSet == 2
                    for iSession = 2
                        sc = nan;

                        if RAND == true
                            randCell = clusterNums{randi(numClusters,1,1)};
                            qCNO = extract.cols(currentData,labels,'quality','cell num',randCell,'session',sessions{iSession});
                            mapCNO = extract.cols(currentData,labels,'rate map','cell num',randCell,'session',sessions{iSession});
                            peakRateCNO = extract.cols(currentData,labels,'peak rate','cell num',randCell,'session',sessions{iSession});
                        else
                            qCNO = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{iSession});
                            mapCNO = extract.cols(currentData,labels,'rate map','cell num',clusterNums{iCluster},'session',sessions{iSession});
                            peakRateCNO = extract.cols(currentData,labels,'peak rate','cell num',clusterNums{iCluster},'session',sessions{iSession});
                        end
                         
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
                            
                            % calculate distance b/t current & previous peak
                            distToCNO = pdist([prevPeakX,prevPeakY; CNOprimaryPeakX,CNOprimaryPeakY]);
                            distStore(iCluster,1) = distToCNO;
                            sc = analyses.spatialCrossCorrelation(mapBL,mapCNO);
                            
                            % for cells that remap
                            if distToCNO > 10 % 
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
                                spatialCorr(iCluster,iSession) = analyses.spatialCrossCorrelation(mapBL,mapCNO);
                                
                                % peak firing in predicted field
                                minDist = nanmin(predictionError(iCluster,:));
                                idx = find(predictionError(iCluster,:) == minDist);
                                if ~isempty(idx)
                                    ratePredicted(iCluster,2) = ratesSorted(idx(1)+1);
                                end
                                
                                % plot CNO map & location of primary peaks
                                if PLOT == true
                                    subplot(1,2,iSession)
                                    hold on
                                    colorMapBRK(mapCNO);
                                    plot(CNOprimaryPeakX,CNOprimaryPeakY,'+','markeredgecolor','k','markersize',15,'linewidth',2.5)
                                    text(5,35,sprintf('Distance = %.2f',distToCNO));
                                    text(5,40,sprintf('Error = %.2f',nanmin(predictionError(iCluster,:))));
                                end
                            end
                        end
                    end
                    
                elseif iDataSet == 3 || iDataSet == 4
                    for iSession = 1:length(sessions)
                        sc = nan;
                        
                        if RAND == true
                            randCell = clusterNums{randi(numClusters,1,1)};
                            qCNO = extract.cols(currentData,labels,'quality','cell num',randCell,'session',sessions{iSession});
                        else
                            qCNO = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{iSession});
                        end
                        
                        if ~strcmpi(qCNO,badQ)
                            numBlocks = 2;
                            
                            if RAND == true
                                posAve = extract.cols(currentData,labels,'p','cell num',clusterNums{randCell},'session',sessions{2});
                                blockLength = floor(size(posAve(:,1),1)/numBlocks);
                                spikes = extract.cols(currentData,labels,'s','cell num',clusterNums{randCell},'session',sessions{2});
                                [s1 t1 x1 y1] = sessionBlocks(1,numBlocks,blockLength,posAve,spikes);
                                mapCNO = analyses.map([t1 x1 y1],s1,'smooth',hippoGlobe.smoothing,'binWidth',hippoGlobe.binWidth,'limits',hippoGlobe.mapLimits);
                                peakRateCNO = mapCNO.peakRate;
                            else
                                posAve = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'session',sessions{2});
                                blockLength = floor(size(posAve(:,1),1)/numBlocks);
                                spikes = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'session',sessions{2});
                                [s1 t1 x1 y1] = sessionBlocks(1,numBlocks,blockLength,posAve,spikes);
                                mapCNO = analyses.map([t1 x1 y1],s1,'smooth',hippoGlobe.smoothing,'binWidth',hippoGlobe.binWidth,'limits',hippoGlobe.mapLimits);
                                peakRateCNO = mapCNO.peakRate;
                            end
                            
                            [~,CNOfields] = analyses.placefield(mapCNO.z,'threshold',0.4,'binWidth',hippoGlobe.binWidth,'minBins',10,'minPeak',0.5);
                            
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
                                sc = analyses.spatialCrossCorrelation(mapBL,mapCNO.z);
                                
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
                                    spatialCorr(iCluster,iSession) = analyses.spatialCrossCorrelation(mapBL,mapCNO.z);
                                    
                                    % peak firing in predicted field
                                    minDist = nanmin(predictionError(iCluster,:));
                                    idx = find(predictionError(iCluster,:) == minDist);
                                    if ~isempty(idx)
                                        ratePredicted(iCluster,2) = ratesSorted(idx(1)+1); 
                                    end                             
         
                                    % plot CNO map & location of primary peaks
                                    if PLOT == true
                                        subplot(1,2,iSession)
                                        hold on
                                        colorMapBRK(mapCNO.z);
                                        plot(CNOprimaryPeakX,CNOprimaryPeakY,'+','markeredgecolor','k','markersize',15,'linewidth',2.5)
                                        text(5,35,sprintf('Distance = %.2f',distToCNO));
                                        text(5,40,sprintf('Error = %.2f',nanmin(predictionError(iCluster,:))));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    if iDataSet == 1
        predictionsAldis = predictionError;
        errorValsAldis = errorVals;
        spatialCorrAldis = spatialCorr;
        ratePredictedAldis = ratePredicted;
        distStoreAldis = distStore;
    elseif iDataSet == 2
        predictionsJas = predictionError;
        errorValsJas = errorVals;
        spatialCorrJas = spatialCorr;
        ratePredictedJas = ratePredicted;
        distStoreJas = distStore;
    elseif iDataSet == 3
        predictionsCont = predictionError;
        errorValsCont = errorVals;
        spatialCorrCont = spatialCorr;
        ratePredictedCont = ratePredicted;
        distStoreCont = distStore;      
    elseif iDataSet == 4
        predictionsBreak = predictionError;
        errorValsBreak = errorVals;
        spatialCorrBreak = spatialCorr;
        ratePredictedBreak = ratePredicted; 
        distStoreBreak = distStore;
    end
end

% pool groups
predictionsTotal = [predictionsAldis; predictionsJas; predictionsCont; predictionsBreak];
errorValsTotal = [errorValsAldis; errorValsJas; errorValsCont; errorValsBreak];
spatialCorrTotal = [spatialCorrAldis; spatialCorrJas; spatialCorrCont; spatialCorrBreak];
ratePredictedTotal = [ratePredictedAldis; ratePredictedJas; ratePredictedCont; ratePredictedBreak];
distStoreTotal = [distStoreAldis; distStoreJas; distStoreCont; distStoreBreak];

%%
peakRate = ratePredictedTotal(:,2);
spatialCorr = spatialCorrTotal(:,2);

err = errorValsTotal(:,2) * hippoGlobe.binWidth;
nanmedian(err)

%%
binRange = 0:2:60;
counts = histcounts(err,binRange) / sum(~isnan(err));
figure;
bar(1:2:59,counts,'barwidth',1,'edgecolor','k','facecolor',[0.5 0.5 0.5])
set(gca,'xtickmode','manual','xtick',0:10:60,'box','off','fontsize',14,'fontweight','bold')
xlabel('Prediction error (cm)')
ylabel('Proportion of neurons')
ylim([0 0.16])

nanmedian(err)
sum(~isnan(err))
