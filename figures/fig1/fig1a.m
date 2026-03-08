load fig1a

figure;
cnt = 1;
for i = 1:3
    for j = 1:4
        subplot(3,4,cnt)
        colorMapBRK(mapStore{i,j});
        text(0,size(mapStore{i,j},1)+3,sprintf('%d Hz',round(peakRates(i,j))))
        cnt = cnt + 1;
    end
end

%%
startup

fieldThresh = 0.2; 
fieldPeak = 0.5; 

rateStore = nan(size(mapStore,1),10,4);

for i = 1:size(mapStore,1)
    map = mapStore{i,1}; % BL1
    [fieldMap, fieldStruct] = analyses.placefield(map,'threshold',fieldThresh,'binWidth',hippoGlobe.binWidth,'minBins',10,'minPeak',fieldPeak);
    
    if length(fieldStruct) > 1
        c = nchoosek(1:length(fieldStruct),2);
        toRemove = zeros(length(fieldStruct),1);
        ovrlap = zeros(size(c,1),1);
        for iField = 1:size(c,1)
            ovrlap(iField) = rectint(fieldStruct(c(iField,1)).bbox,fieldStruct(c(iField,2)).bbox + [-1 -1 2 2]);
            if ovrlap(iField) > 1
                a1 = fieldStruct(c(iField,1)).area;
                a2 = fieldStruct(c(iField,2)).area;
                if a1 > a2
                    toRemove(c(iField,2)) = 1;
                elseif a2 > a1
                    toRemove(c(iField,1)) = 1;
                end
            end
        end
        fieldStruct = fieldStruct(~toRemove);
        %             fsMapNew = fieldMap;
        remList = find(toRemove);
        for iField = 1:length(remList)
            fieldMap(fieldMap == remList(iField)) = 0;
        end
    end
    
    if ~isempty(fieldStruct)
        numFields = length(fieldStruct);
        if numFields > 1
            rates = []; pixList = {};
            % get peak rates and pixel lists for each field
            for iField = 1:numFields
                rates(iField) = fieldStruct(1,iField).peak;
                pixList{iField} = fieldStruct(iField).PixelIdxList;
            end
            % sort by peak rate and preserve order for other sessions
            [~,sortInd] = sort(rates,'descend');
            % store peak rates
            for iField = 1:numFields
                rateStore(i,iField,1) = rates(sortInd(iField));
            end
            
            %%
            map = mapStore{i,2}; % CNO1        
            sortedRates = [];
            % for each BL field, store peak rate of those pixels in CNO session
            for iField = 1:numFields
                % if more than half of these pixels are nan then forget it
                if sum(isnan(map(pixList{sortInd(iField)}))) > (0.5*numel(map(pixList{sortInd(iField)})))
                    sortedRates(iField) = nan;
                else
                    sortedRates(iField) = max(nanmax(map(pixList{sortInd(iField)})));
                end
            end
            for iField = 1:numFields
                rateStore(i,iField,2) = sortedRates(iField);
            end
        end
        
        %%
        map = mapStore{i,3}; % BL2
        sortedRates = [];
        % for each BL field, store peak rate of those pixels in CNO session
        for iField = 1:numFields
            % if more than half of these pixels are nan then forget it
            if sum(isnan(map(pixList{sortInd(iField)}))) > (0.5*numel(map(pixList{sortInd(iField)})))
                sortedRates(iField) = nan;
            else
                sortedRates(iField) = max(nanmax(map(pixList{sortInd(iField)})));
            end
        end
        for iField = 1:numFields
            rateStore(i,iField,3) = sortedRates(iField);
        end
        
        %%
        map = mapStore{i,4}; % CNO2
        sortedRates = [];
        % for each BL field, store peak rate of those pixels in CNO session
        for iField = 1:numFields
            % if more than half of these pixels are nan then forget it
            if sum(isnan(map(pixList{sortInd(iField)}))) > (0.5*numel(map(pixList{sortInd(iField)})))
                sortedRates(iField) = nan;
            else
                sortedRates(iField) = max(nanmax(map(pixList{sortInd(iField)})));
            end
        end
        for iField = 1:numFields
            rateStore(i,iField,4) = sortedRates(iField);
        end
    end
end

%%
figure;
cnt = 1;
for i = 1:size(rateStore,1)
    tempRates = rateStore(i,:,:);
    tempRates = reshape(tempRates, 10, 4).';
    
    inds = find(sum(isnan(tempRates)) == 4);
    n = inds(1);
    
    for j = 1:4
        subplot(3,4,cnt)
        h = bar(1:n,tempRates(j,1:n));
        set(h,'barwidth',0.8,'facecolor',[.5 .5 .5])
       ylim([0 round(nanmax(nanmax(tempRates)))])
        set(gca,'xtickmode','manual','box','off','fontsize',8)
        xlabel('Subfield')
        ylabel('Peak rate (Hz)')
        cnt = cnt + 1;
    end
end