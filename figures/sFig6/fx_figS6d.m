
dataStore = [];
load fx_figS6d_indRealign
load fx_figS6d_subfieldRate

PLOT = false; % turn on to plot excitation maps and corresponding place cell rate maps
RAND = false;
    
%%
for iSim = 1:numSims

    iBL_CA3 = iBL_CA3_all{iSim,1};
    pBL_CA3 = pBL_CA3_all{iSim,1};
    iCNO_CA3 = iCNO_CA3_all{iSim,1};
    pCNO_CA3 = pCNO_CA3_all{iSim,1};
    CA3gridInputs = CA3gridInputs_all{iSim,1};
    CA3weights = CA3weights_all{iSim,1};
    
    % stores cell num, spatial corr CA3 vs. MEC, prediction from MEC, sim num
    simStore = nan(200,5);
    
    for iPlace = 1:200 
        
        predictionError = nan(1,3);
        simStore(iPlace,1,:) = iPlace;
        simStore(iPlace,2,:) = iSim;
        
        %% BL
        % total MEC --> CA3
        tmp = zeros(100);
        for iGridInput = 1:size(CA3gridInputs.randGridStore,2)
            tmp = tmp + ...
                MEC_BL.GridPool(:,:,CA3gridInputs.randGridStore(iPlace,iGridInput)) * ...
                CA3weights.WeightMatrix(iGridInput,iPlace);
        end
        BL.MECexcitation = tmp;
        MECmapBL = BL.MECexcitation;
        
        % prediction on MEC map
        [~,fieldsMEC] = analyses.placefield(MECmapBL,'threshold',0.1,'minPeak',0.5);
        
        if ~isempty(fieldsMEC)
            peaks = [];
            for iFieldNum = 1:size(fieldsMEC,2)
                peaks(iFieldNum) = fieldsMEC(1,iFieldNum).peak;
            end
            primary = find(peaks == max(peaks));
            if length(primary) > 1
                primary = primary(1);
            end
            MECprimaryPeakX = fieldsMEC(1,primary).peakX;
            MECprimaryPeakY = fieldsMEC(1,primary).peakY;
            
            map = MECmapBL;
            map(isnan(map)) = 0;
            regmax = imregionalmax(map);
            [y x] = find(regmax == 1);
            
            meanRates = [];
            for j = 1:length(x)
                meanRates(j) = MECmapBL(y(j),x(j));
            end
            
            % find location of second & third peaks in BL & plot
            ratesSorted = sort(meanRates,'descend');
            if length(ratesSorted) > 1
                secondIdx = find(ratesSorted(2) == meanRates);
            end
            if length(ratesSorted) > 2
                thirdIdx = find(ratesSorted(3) == meanRates);
            end
            if length(ratesSorted) > 3
                fourthIdx = find(ratesSorted(4) == meanRates);
            end
            
            if PLOT == true
                figure('name',sprintf('Cell %d',iPlace))
                subplot(221)
                hold on
                colorMapBRK(MECmapBL,'bar','on','clrmap','jet');
                plot(MECprimaryPeakX,MECprimaryPeakY,'+','markeredgecolor','k','markersize',15,'linewidth',2.5)
                % title(sprintf('SC = %.2f',scMEC))
                
                if length(ratesSorted) > 1
                    plot(x(secondIdx),y(secondIdx),'*','markeredgecolor','m','markersize',15,'linewidth',2.5)
                end
                if length(ratesSorted) > 2
                    plot(x(thirdIdx),y(thirdIdx),'x','markeredgecolor','b','markersize',15,'linewidth',2.5)
                end
                if length(ratesSorted) > 3
                    plot(x(fourthIdx),y(fourthIdx),'d','markeredgecolor','y','markersize',15,'linewidth',2.5)
                end
            end
            
            %% CNO
            % total MEC --> CA3
            if RAND == true
                iPlaceRand = randi(200,1,1);
                tmp = zeros(100);
                for iGridInput = 1:size(CA3gridInputs.randGridStore,2)
                    tmp = tmp + ...
                        MEC_CNO.GridPool(:,:,CA3gridInputs.randGridStore(iPlaceRand,iGridInput)) * ...
                        CA3weights.WeightMatrix(iGridInput,iPlaceRand);
                end
                CNO.MECexcitation = tmp;
                MECmapCNO = CNO.MECexcitation;
            else
                tmp = zeros(100);
                for iGridInput = 1:size(CA3gridInputs.randGridStore,2)
                    tmp = tmp + ...
                        MEC_CNO.GridPool(:,:,CA3gridInputs.randGridStore(iPlace,iGridInput)) * ...
                        CA3weights.WeightMatrix(iGridInput,iPlace);
                end
                CNO.MECexcitation = tmp;
                MECmapCNO = CNO.MECexcitation;
            end
            
            % prediction on MEC map
            [~,fieldsMEC] = analyses.placefield(MECmapCNO,'threshold',0.1,'minPeak',0.5);
            
            if ~isempty(fieldsMEC)
                peaks = [];
                for iFieldNum = 1:size(fieldsMEC,2)
                    peaks(iFieldNum) = fieldsMEC(1,iFieldNum).peak;
                end
                primary = find(peaks == max(peaks));
                if length(primary) > 1
                    primary = primary(1);
                end
                MECprimaryPeakX_CNO = fieldsMEC(1,primary).peakX;
                MECprimaryPeakY_CNO = fieldsMEC(1,primary).peakY;
                
                distToBL = pdist([MECprimaryPeakX,MECprimaryPeakY; MECprimaryPeakX_CNO,MECprimaryPeakY_CNO]);
                scMEC = analyses.spatialCrossCorrelation(MECmapBL,MECmapCNO);
                
                % find location of second & third peaks in CNO & plot
                ratesSorted = sort(meanRates,'descend');
                if length(ratesSorted) > 1
                    secondIdx = find(ratesSorted(2) == meanRates);
                    predictionError(1) = pdist([x(secondIdx),y(secondIdx); MECprimaryPeakX_CNO,MECprimaryPeakY_CNO]);
                end
                if length(ratesSorted) > 2
                    thirdIdx = find(ratesSorted(3) == meanRates);
                    predictionError(2) = pdist([x(thirdIdx),y(thirdIdx); MECprimaryPeakX_CNO,MECprimaryPeakY_CNO]);
                end
                if length(ratesSorted) > 3
                    fourthIdx = find(ratesSorted(4) == meanRates);
                    predictionError(3) = pdist([x(fourthIdx),y(fourthIdx); MECprimaryPeakX_CNO,MECprimaryPeakY_CNO]);
                end
                
                % minimum prediction for each area
                minMEC = nanmin(predictionError(1:3));
                
                % store values
                simStore(iPlace,3) = scMEC;
                simStore(iPlace,4) = distToBL;
                simStore(iPlace,5) = minMEC;
                
                if PLOT == true
                    subplot(222)
                    hold on
                    colorMapBRK(MECmapCNO,'bar','on','clrmap','jet');
                    plot(MECprimaryPeakX_CNO,MECprimaryPeakY_CNO,'+','markeredgecolor','k','markersize',15,'linewidth',2.5)
                    title(sprintf('SC = %.2f',scMEC))
                    text(0,120,sprintf('Dist to BL = %.2f',distToBL),'fontsize',14)
                    text(0,140,sprintf('MEC dist = %.2f',minMEC),'fontsize',14)
                    
                    mapBL = general.smooth(pBL_CA3(:,:,iPlace),6.67);
                    mapCNO = general.smooth(pCNO_CA3(:,:,iPlace),6.67);
                   
                    subplot(223)      
                    colorMapBRK(mapBL);
                    
                    subplot(224)
                    colorMapBRK(mapCNO);
                end
            end
        end
    end
    dataStore = [dataStore; simStore];
end

%
nanmedian(dataStore)
sum(~isnan(dataStore))
dataStoreShift = dataStore;

%% CDF
figure;
hold on

[f,x,flo,fup] = ecdf(dataStoreShift(:,3));
shadedplot(x',flo',fup','w','k');
C = cdfplot(dataStoreShift(:,3));
set(C,'color','k','linewidth',0.5)

title ''
grid off
box off
ylabel 'Proportion'
xlabel 'Spatial correlation'
set(gca,'fontsize',14,'fontw','bold')
xlim([-0.6 1])
