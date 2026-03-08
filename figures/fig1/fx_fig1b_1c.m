function fx_fig1b_1c

global hippoGlobe

% check inputs and load data
cellType = 'grid'; 
PLOT = false;
RAND = false;

for iDataSet = 1:3
    if iDataSet == 1
        load masterMat_MEC1
        RECT
    elseif iDataSet == 2
        load masterMat_MEC2
        BIGBOX
    elseif iDataSet == 3
        load masterMat_MEC3
        BOX
    end
    
    raw = dataOutput;
    
    %% store names of groups and sessions
    sessions = unique(extract.cols(raw,labels,'session'),'stable');
    numSesh = numel(sessions);
    
    badQ = '3';
    offQ = '4';
    
    if iDataSet == 1
        % switch numbers to strings
        for iRow = 1:size(raw,1)
            raw{iRow,strcmpi('exptHalf',labels)} = num2str(raw{iRow,strcmpi('exptHalf',labels)});
        end
        
        expHalves = unique(extract.cols(raw,labels,'exptHalf'),'sorted');
        numExpHalves = numel(expHalves);
        raw = cleanUpQualityMEC(raw,labels,expHalves,numExpHalves,sessions,numSesh,badQ,offQ);
    elseif iDataSet == 2
        raw = extract.cleanUpQuality(raw,labels,sessions,numSesh,badQ,offQ);
    elseif iDataSet == 3
        raw = extract.cleanUpQuality(raw,labels,sessions,numSesh,badQ,offQ);
    end
    
    %% functional thresholds
    fieldThresh = 0.2; 
    fieldPeak = 0.5; 
    
    load MECshuffled
    thresh.FR = 10;
    thresh.HD = shuffled.HD;
    thresh.grid = shuffled.grid;
    thresh.border = shuffled.border;
    thresh.spatial = shuffled.spatInfo;
    thresh.theta = 4;
    
    %% filter for functional cellType
    if iDataSet == 1
        if exist('cellType','var')
            if sum(strcmpi(cellType,{'excitatory','HD','grid','border','spatial'}))
                toRemove = extract.cols(raw,labels,'cell num','exptHalf',expHalves{1},'session',sessions{1},'mean rate','>=',thresh.FR);
                raw = extract.rows(raw,labels,'remove','cell num',toRemove);
            end
            if strcmpi(cellType,'inhibitory')
                toRemove = extract.cols(raw,labels,'cell num','exptHalf',expHalves{1},'session',sessions{1},'mean rate','<',thresh.FR);
                raw = extract.rows(raw,labels,'remove','cell num',toRemove);
            elseif strcmpi(cellType,'HD')
                toRemove = raw(strcmpi(raw(:,strcmpi('session',labels)),sessions{1}) ...
                    & strcmpi(raw(:,strcmpi('exptHalf',labels)),expHalves{1}) ...
                    & (cell2mat(raw(:,strcmpi('mean vector length',labels))) < thresh.HD ...
                    | cellfun(@isnan,raw(:,strcmpi('mean vector length',labels)))), ...
                    strcmpi('cell num',labels));
                raw = extract.rows(raw,labels,'remove','cell num',toRemove);
            elseif strcmpi(cellType,'grid')
                toRemove = raw(strcmpi(raw(:,strcmpi('session',labels)),sessions{1}) ...
                    & strcmpi(raw(:,strcmpi('exptHalf',labels)),expHalves{1}) ...
                    & (cell2mat(raw(:,strcmpi('grid score',labels))) < thresh.grid ...
                    | cellfun(@isnan,raw(:,strcmpi('grid score',labels)))), ...
                    strcmpi('cell num',labels));
                raw = extract.rows(raw,labels,'remove','cell num',toRemove);
            elseif strcmpi(cellType,'border')
                toRemove = raw(strcmpi(raw(:,strcmpi('session',labels)),sessions{1}) ...
                    & strcmpi(raw(:,strcmpi('exptHalf',labels)),expHalves{1}) ...
                    & (cell2mat(raw(:,strcmpi('border score',labels))) < thresh.border ...
                    | cellfun(@isnan,raw(:,strcmpi('border score',labels)))), ...
                    strcmpi('cell num',labels));
                raw = extract.rows(raw,labels,'remove','cell num',toRemove);
                % high spatial info
                toRemove = raw(strcmpi(raw(:,strcmpi('session',labels)),sessions{1}) ...
                    & strcmpi(raw(:,strcmpi('exptHalf',labels)),expHalves{1}) ...
                    & cell2mat(raw(:,strcmpi('spatial info',labels))) < thresh.spatial, ...
                    strcmpi('cell num',labels));
                raw = extract.rows(raw,labels,'remove','cell num',toRemove);
            elseif strcmpi(cellType,'spatial')
                toRemove = raw(strcmpi(raw(:,strcmpi('session',labels)),sessions{1}) ...
                    & strcmpi(raw(:,strcmpi('exptHalf',labels)),expHalves{1}) ...
                    & (cell2mat(raw(:,strcmpi('spatial info',labels))) < thresh.spatial ...
                    | cell2mat(raw(:,strcmpi('grid score',labels))) >= thresh.grid ...
                    | cell2mat(raw(:,strcmpi('border score',labels))) >= thresh.border), ...
                    strcmpi('cell num',labels));
                raw = extract.rows(raw,labels,'remove','cell num',toRemove);
            end
        end
        
    elseif iDataSet == 2 || iDataSet == 3
        if exist('cellType','var')
            if sum(strcmpi(cellType,{'excitatory','HD','grid','border','spatial'}))
                toRemove = extract.cols(raw,labels,'cell num','session',sessions{1},'mean rate','>=',thresh.FR);
                raw = extract.rows(raw,labels,'remove','cell num',toRemove);
            end
            if strcmpi(cellType,'inhibitory')
                toRemove = extract.cols(raw,labels,'cell num','session',sessions{1},'mean rate','<',thresh.FR);
                raw = extract.rows(raw,labels,'remove','cell num',toRemove);
            elseif strcmpi(cellType,'HD')
                toRemove = raw(strcmpi(raw(:,strcmpi('session',labels)),sessions{1}) ...
                    & (cell2mat(raw(:,strcmpi('mean vector length',labels))) < thresh.HD ...
                    | cellfun(@isnan,raw(:,strcmpi('mean vector length',labels)))), ...
                    strcmpi('cell num',labels));
                raw = extract.rows(raw,labels,'remove','cell num',toRemove);
            elseif strcmpi(cellType,'grid')
                toRemove = raw(strcmpi(raw(:,strcmpi('session',labels)),sessions{1}) ...
                    & (cell2mat(raw(:,strcmpi('grid score',labels))) < thresh.grid ...
                    | cellfun(@isnan,raw(:,strcmpi('grid score',labels)))), ...
                    strcmpi('cell num',labels));
                raw = extract.rows(raw,labels,'remove','cell num',toRemove);
            elseif strcmpi(cellType,'border')
                toRemove = raw(strcmpi(raw(:,strcmpi('session',labels)),sessions{1}) ...
                    & (cell2mat(raw(:,strcmpi('border score',labels))) < thresh.border ...
                    | cellfun(@isnan,raw(:,strcmpi('border score',labels)))), ...
                    strcmpi('cell num',labels));
                raw = extract.rows(raw,labels,'remove','cell num',toRemove);
                % high spatial info
                toRemove = raw(strcmpi(raw(:,strcmpi('session',labels)),sessions{1}) ...
                    & cell2mat(raw(:,strcmpi('spatial info',labels))) < thresh.spatial, ...
                    strcmpi('cell num',labels));
                raw = extract.rows(raw,labels,'remove','cell num',toRemove);
            elseif strcmpi(cellType,'spatial')
                toRemove = raw(strcmpi(raw(:,strcmpi('session',labels)),sessions{1}) ...
                    & (cell2mat(raw(:,strcmpi('spatial info',labels))) < thresh.spatial ...
                    | cell2mat(raw(:,strcmpi('grid score',labels))) >= thresh.grid ...
                    | cell2mat(raw(:,strcmpi('border score',labels))) >= thresh.border), ...
                    strcmpi('cell num',labels));
                raw = extract.rows(raw,labels,'remove','cell num',toRemove);
            end
        end
    end
    
    %% separate functional groups
    if iDataSet == 1
        groupNames = unique(extract.cols(raw,labels,'group'),'sorted');
        numGroups = numel(groupNames);
        Con = raw(strcmpi(raw(:,strcmpi('group',labels)),'con'),:);
        Exp = raw(strcmpi(raw(:,strcmpi('group',labels)),'exp'),:);
        
        % con
        mapStore = nan(23,30,500,4,1); % pix x pix x num cells x session x group
        % exp
        mapStore = nan(23,30,500,4,2);
        
    elseif iDataSet == 2
        groupNames = unique(extract.cols(raw,labels,'group'),'sorted');
        numGroups = numel(groupNames);
        Con = raw(strcmpi(raw(:,strcmpi('group',labels)),'con'),:);
        
        % con
        mapStore = nan(25,25,500,5,1); % pix x pix x num cells x session x group
        
    elseif iDataSet == 3
        groupNames = {'Con','Exp'};
        numGroups = numel(groupNames);
        Con = extract.cols(raw,labels,':','genotype','exp','injection','sal');
        Exp = extract.cols(raw,labels,':','genotype','exp','injection','cno');
        
        % con
        mapStore = nan(30,30,500,5,1); % pix x pix x num cells x session x group
        % exp
        mapStore = nan(30,30,500,5,2);
    end
    
    %% initialize
    rateStore = nan(50,10,4,4);            % cell x # grid fields x session x group
    
    if iDataSet == 1
        for iGroup = 1:numGroups
            eval(sprintf('currentData = %s;',groupNames{iGroup}))
            clusterNums = unique(extract.cols(currentData,labels,'cell num'),'stable');
            numClusters = length(clusterNums);
            
            for iCluster = 1:numClusters
                
                % quality for each piece
                q1 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{1});
                q3 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{3});
                q4 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{4});
                
                % BL1
                % get rate map and find fields
                if ~strcmpi(q1,badQ)
                    map = currentData{strcmpi(currentData(:,strcmpi('cell num',labels)),clusterNums{iCluster}) ...
                        & strcmpi(currentData(:,strcmpi('exptHalf',labels)),expHalves{1}) ...
                        & strcmpi(currentData(:,strcmpi('session',labels)),sessions{1}), ...
                        strcmpi('rate map',labels)};
                    mapStore(:,:,iCluster,1,iGroup) = map;
                    [fieldMap, fieldStruct] = analyses.placefield(map,'threshold',fieldThresh,'binWidth',hippoGlobe.binWidth,'minBins',2,'minPeak',fieldPeak);
                    
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
                                rateStore(iCluster,iField,1,iGroup) = rates(sortInd(iField));
                            end
                            
                            if PLOT == true
                                figure;
                                subplot(351)
                                colorMapBRK(map);
                                title(sprintf('Cell ID = %s',clusterNums{iCluster}))
                                hold on
                                for iField = 1:length(fieldStruct)
                                    rectangle('position',fieldStruct(iField).bbox,'linew',3)
                                end
                                subplot(356)
                                bar(1:length(rates),rateStore(iCluster,1:length(rates),1,iGroup))
                            end

                            % CNO1
                            % C1
                            if ~strcmpi(q3,badQ)
                                p2 = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{3});
                                s2 = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{3});
                            end
                            
                            % D1
                            if ~strcmpi(q4,badQ)
                                p3 = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{4});
                                s3 = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{4});
                            end
                            
                            % merge sessions
                            if ~strcmpi(q3,badQ) && ~strcmpi(q4,badQ)
                                offset = mean(diff(p2(:,1))) - min(p3(:,1)) + p2(end,1);
                                sCNOday1 = [s2; s3+offset];
                                pCNOday1 = [p2; p3+offset];
                                  
                                map = analyses.map(pCNOday1,sCNOday1,'smooth',hippoGlobe.smoothing,'binWidth',hippoGlobe.binWidth,'limits',hippoGlobe.mapLimits);
                                map = map.z;
                                mapStore(:,:,iCluster,2,iGroup) = map;
                                
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
                                    rateStore(iCluster,iField,2,iGroup) = sortedRates(iField);
                                end
                                
                                if PLOT == true
                                    subplot(352)
                                    colorMapBRK(map);
                                    hold on
                                    for iField = 1:length(fieldStruct)
                                        rectangle('position',fieldStruct(iField).bbox,'linew',3)
                                    end
                                    subplot(357)
                                    bar(1:length(rates),rateStore(iCluster,1:length(rates),2,iGroup));
                                end
                                
                                q6 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{1});
                                if RAND == true
                                    randCell = randi(numClusters);
                                    q8 = extract.cols(currentData,labels,'quality','cell num',clusterNums{randCell},'exptHalf',expHalves{2},'session',sessions{3});
                                    q9 = extract.cols(currentData,labels,'quality','cell num',clusterNums{randCell},'exptHalf',expHalves{2},'session',sessions{4});
                                else
                                    q8 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{3});
                                    q9 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{4});
                                end
                            end
                            
                            % BL2
                            if ~strcmpi(q6,badQ)
                                map = currentData{strcmpi(currentData(:,strcmpi('cell num',labels)),clusterNums{iCluster}) ...
                                    & strcmpi(currentData(:,strcmpi('exptHalf',labels)),expHalves{2}) ...
                                    & strcmpi(currentData(:,strcmpi('session',labels)),sessions{1}), ...
                                    strcmpi('rate map',labels)};
                                mapStore(:,:,iCluster,3,iGroup) = map;
                                
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
                                    rateStore(iCluster,iField,3,iGroup) = sortedRates(iField);
                                end
                                
                                if PLOT == true
                                    subplot(353)
                                    colorMapBRK(map);
                                    hold on
                                    for iField = 1:length(fieldStruct)
                                        rectangle('position',fieldStruct(iField).bbox,'linew',3)
                                    end
                                    subplot(358)
                                    bar(1:length(rates),rateStore(iCluster,1:length(rates),3,iGroup));
                                end
                            end
                            
                            % C2
                            if ~strcmpi(q8,badQ)
                                if RAND == true
                                    p5 = extract.cols(currentData,labels,'p','cell num',clusterNums{randCell},'exptHalf',expHalves{2},'session',sessions{3});
                                    s5 = extract.cols(currentData,labels,'s','cell num',clusterNums{randCell},'exptHalf',expHalves{2},'session',sessions{3});
                                else
                                    p5 = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{3});
                                    s5 = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{3});
                                end
                            end
                            
                            % D2
                            if ~strcmpi(q9,badQ)
                                if RAND == true
                                    p6 = extract.cols(currentData,labels,'p','cell num',clusterNums{randCell},'exptHalf',expHalves{2},'session',sessions{4});
                                    s6 = extract.cols(currentData,labels,'s','cell num',clusterNums{randCell},'exptHalf',expHalves{2},'session',sessions{4});
                                else
                                    p6 = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{4});
                                    s6 = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{4});
                                end
                            end
                            
                            if ~strcmpi(q8,badQ) && ~strcmpi(q9,badQ)
                                offset = mean(diff(p5(:,1))) - min(p6(:,1)) + p5(end,1);
                                sCNOday2 = [s5; s6+offset];
                                pCNOday2 = [p5; p6+offset];
                                
                                map = analyses.map(pCNOday2,sCNOday2,'smooth',hippoGlobe.smoothing,'binWidth',hippoGlobe.binWidth,'limits',hippoGlobe.mapLimits);
                                map = map.z;
                                mapStore(:,:,iCluster,4,iGroup) = map;
                                
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
                                    rateStore(iCluster,iField,4,iGroup) = sortedRates(iField);
                                end
                                
                                if PLOT == true
                                    subplot(354)
                                    colorMapBRK(map);
                                    hold on
                                    for iField = 1:length(fieldStruct)
                                        rectangle('position',fieldStruct(iField).bbox,'linew',3)
                                    end
                                    subplot(359)
                                    bar(1:length(rates),rateStore(iCluster,1:length(rates),4,iGroup));
                                end
                            end
                        end
                    end
                end
            end
        end
        
        
    elseif iDataSet == 2
    
        for iGroup = 1:numGroups
            currentData = Con;
            clusterNums = unique(extract.cols(currentData,labels,'cell num'),'stable');
            numClusters = length(clusterNums);
            
            for iCluster = 1:numClusters
                
                % quality for each piece
                q1 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{1});
                q2 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{2});
                q3 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{3});
                
                if RAND == true
                    randCell = randi(numClusters);
                    q4 = extract.cols(currentData,labels,'quality','cell num',clusterNums{randCell},'session',sessions{4});
                else
                    q4 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{4});
                end
                
                % BL1
                % get rate map and find fields
                if ~strcmpi(q1,badQ)
                    map = currentData{strcmpi(currentData(:,strcmpi('cell num',labels)),clusterNums{iCluster}) ...
                        & strcmpi(currentData(:,strcmpi('session',labels)),sessions{1}), ...
                        strcmpi('rate map',labels)};
                    mapStore(:,:,iCluster,1,iGroup) = map;
                    [fieldMap, fieldStruct] = analyses.placefield(map,'threshold',fieldThresh,'binWidth',hippoGlobe.binWidth,'minBins',2,'minPeak',fieldPeak);
                    
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
                                rateStore(iCluster,iField,1,iGroup) = rates(sortInd(iField));
                            end
                            
                            if PLOT == true
                                figure;
                                subplot(351)
                                colorMapBRK(map);
                                title(sprintf('Cell ID = %s',clusterNums{iCluster}))
                                hold on
                                for iField = 1:length(fieldStruct)
                                    rectangle('position',fieldStruct(iField).bbox,'linew',3)
                                end
                                subplot(356)
                                bar(1:length(rates),rateStore(iCluster,1:length(rates),1,iGroup))
                            end
                            
                            numBlocks = 4;
                            % CNO1
                            if ~strcmpi(q2,badQ)
                                posAve = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'session',sessions{2});
                                blockLength = floor(size(posAve(:,1),1)/numBlocks);
                                spikes = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'session',sessions{2});
                                [s2 t2 x2 y2] = sessionBlocks(2,numBlocks,blockLength,posAve,spikes);
                                [s3 t3 x3 y3] = sessionBlocks(3,numBlocks,blockLength,posAve,spikes);
                                
                                % merge the sessions
                                offset = mean(diff(t2)) - min(t3) + t2(end);
                                sCNOday1 = [s2; s3+offset];
                                pCNOday1 = [t2 x2 y2; t3+offset x3 y3];
                                
                                map = analyses.map(pCNOday1,sCNOday1,'smooth',hippoGlobe.smoothing,'binWidth',hippoGlobe.binWidth,'limits',hippoGlobe.mapLimits);
                                map = map.z;
                                mapStore(:,:,iCluster,2,iGroup) = map;
                                
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
                                    rateStore(iCluster,iField,2,iGroup) = sortedRates(iField);
                                end
                                
                                if PLOT == true
                                    subplot(352)
                                    colorMapBRK(map);
                                    hold on
                                    for iField = 1:length(fieldStruct)
                                        rectangle('position',fieldStruct(iField).bbox,'linew',3)
                                    end
                                    subplot(357)
                                    bar(1:length(rates),rateStore(iCluster,1:length(rates),2,iGroup))
                                end
                            end
                            
                            % BL2
                            if ~strcmpi(q3,badQ)
                                map = currentData{strcmpi(currentData(:,strcmpi('cell num',labels)),clusterNums{iCluster}) ...
                                    & strcmpi(currentData(:,strcmpi('session',labels)),sessions{3}), ...
                                    strcmpi('rate map',labels)};
                                mapStore(:,:,iCluster,3,iGroup) = map;
                                
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
                                    rateStore(iCluster,iField,3,iGroup) = sortedRates(iField);
                                end
                                
                                if PLOT == true
                                    subplot(353)
                                    colorMapBRK(map);
                                    hold on
                                    for iField = 1:length(fieldStruct)
                                        rectangle('position',fieldStruct(iField).bbox,'linew',3)
                                    end
                                    subplot(358)
                                    bar(1:length(rates),rateStore(iCluster,1:length(rates),3,iGroup))
                                end
                            end
                            
                            % CNO2
                            if ~strcmpi(q4,badQ)
                                if RAND == true
                                    % break the CNO session
                                    posAve = extract.cols(currentData,labels,'p','cell num',clusterNums{randCell},'session',sessions{4});
                                    blockLength = floor(size(posAve(:,1),1)/numBlocks);
                                    spikes = extract.cols(currentData,labels,'s','cell num',clusterNums{randCell},'session',sessions{4});
                                    [s5 t5 x5 y5] = sessionBlocks(2,numBlocks,blockLength,posAve,spikes);
                                    [s6 t6 x6 y6] = sessionBlocks(3,numBlocks,blockLength,posAve,spikes);
                                else
                                    % break the CNO session
                                    posAve = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'session',sessions{4});
                                    blockLength = floor(size(posAve(:,1),1)/numBlocks);
                                    spikes = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'session',sessions{4});
                                    [s5 t5 x5 y5] = sessionBlocks(2,numBlocks,blockLength,posAve,spikes);
                                    [s6 t6 x6 y6] = sessionBlocks(3,numBlocks,blockLength,posAve,spikes);
                                end
                                
                                % merge the sessions
                                offset = mean(diff(t5)) - min(t6) + t5(end);
                                sCNOday2 = [s5; s6+offset];
                                pCNOday2 = [t5 x5 y5; t6+offset x6 y6];
                                
                                map = analyses.map(pCNOday2,sCNOday2,'smooth',hippoGlobe.smoothing,'binWidth',hippoGlobe.binWidth,'limits',hippoGlobe.mapLimits);
                                map = map.z;
                                mapStore(:,:,iCluster,4,iGroup) = map;
                                
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
                                    rateStore(iCluster,iField,4,iGroup) = sortedRates(iField);
                                end
                                
                                if PLOT == true
                                    subplot(354)
                                    colorMapBRK(map);
                                    hold on
                                    for iField = 1:length(fieldStruct)
                                        rectangle('position',fieldStruct(iField).bbox,'linew',3)
                                    end
                                    subplot(359)
                                    bar(1:length(rates),rateStore(iCluster,1:length(rates),4,iGroup))
                                end
                            end
                        end
                    end
                end
            end
        end
        
    elseif iDataSet == 3
        for iGroup = 1:numGroups
            eval(sprintf('currentData = %s;',groupNames{iGroup}))
            clusterNums = unique(extract.cols(currentData,labels,'cell num'),'stable');
            numClusters = length(clusterNums);
            
            for iCluster = 1:numClusters
                
                % quality for each piece
                q1 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{1});
                q2 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{2});
                q3 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{3});
                q4 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{4});
                
                if RAND == true
                    randCell = randi(numClusters);
                    q5 = extract.cols(currentData,labels,'quality','cell num',clusterNums{randCell},'session',sessions{5});
                    q6 = extract.cols(currentData,labels,'quality','cell num',clusterNums{randCell},'session',sessions{6});
                else
                    q5 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{5});
                    q6 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{6});
                end
                
                % BL1
                % get rate map and find fields
                if ~strcmpi(q1,badQ)
                    map = currentData{strcmpi(currentData(:,strcmpi('cell num',labels)),clusterNums{iCluster}) ...
                        & strcmpi(currentData(:,strcmpi('session',labels)),sessions{1}), ...
                        strcmpi('rate map',labels)};
                    mapStore(:,:,iCluster,1,iGroup) = map;
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
                                rateStore(iCluster,iField,1,iGroup) = rates(sortInd(iField));
                            end
                            
                            if PLOT == true
                                figure;
                                subplot(351)
                                colorMapBRK(map);
                                title(sprintf('Cell ID = %s',clusterNums{iCluster}))
                                hold on
                                for iField = 1:length(fieldStruct)
                                    rectangle('position',fieldStruct(iField).bbox,'linew',3)
                                end
                                subplot(356)
                                bar(1:length(rates),rateStore(iCluster,1:length(rates),1,iGroup))
                            end
                            
                            % CNO1
                            numBlocks = 2;
                            % break the CNO sessions
                            % 2nd half of CNO1
                            if ~strcmpi(q2,badQ)
                                posAve = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'session',sessions{2});
                                blockLength = floor(size(posAve(:,1),1)/numBlocks);
                                spikes = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'session',sessions{2});
                                [s2 t2 x2 y2] = sessionBlocks(2,numBlocks,blockLength,posAve,spikes);
                            end
                            
                            % 1st half of CNO1
                            if ~strcmpi(q3,badQ)
                                posAve = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'session',sessions{3});
                                blockLength = floor(size(posAve(:,1),1)/numBlocks);
                                spikes = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'session',sessions{3});
                                [s3 t3 x3 y3] = sessionBlocks(1,numBlocks,blockLength,posAve,spikes);
                            end
                            
                            % merge the sessions
                            if ~strcmpi(q2,badQ) && ~strcmpi(q3,badQ)
                                offset = mean(diff(t2)) - min(t3) + t2(end);
                                sCNOday1 = [s2; s3+offset];
                                pCNOday1 = [t2 x2 y2; t3+offset x3 y3];
                                
                                map = analyses.map(pCNOday1,sCNOday1,'smooth',hippoGlobe.smoothing,'binWidth',hippoGlobe.binWidth,'limits',hippoGlobe.mapLimits);
                                map = map.z;
                                mapStore(:,:,iCluster,2,iGroup) = map;
                                
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
                                    rateStore(iCluster,iField,2,iGroup) = sortedRates(iField);
                                end
                                
                                if PLOT == true
                                    subplot(352)
                                    colorMapBRK(map);
                                    hold on
                                    for iField = 1:length(fieldStruct)
                                        rectangle('position',fieldStruct(iField).bbox,'linew',3)
                                    end
                                    subplot(357)
                                    bar(1:length(rates),rateStore(iCluster,1:length(rates),2,iGroup))
                                end
                            end
                            
                            % BL2
                            if ~strcmpi(q3,badQ)
                                map = currentData{strcmpi(currentData(:,strcmpi('cell num',labels)),clusterNums{iCluster}) ...
                                    & strcmpi(currentData(:,strcmpi('session',labels)),sessions{4}), ...
                                    strcmpi('rate map',labels)};
                                mapStore(:,:,iCluster,3,iGroup) = map;
                                
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
                                    rateStore(iCluster,iField,3,iGroup) = sortedRates(iField);
                                end
                                
                                if PLOT == true
                                    subplot(353)
                                    colorMapBRK(map);
                                    hold on
                                    for iField = 1:length(fieldStruct)
                                        rectangle('position',fieldStruct(iField).bbox,'linew',3)
                                    end
                                    subplot(358)
                                    bar(1:length(rates),rateStore(iCluster,1:length(rates),3,iGroup))
                                end
                            end
                            
                            % CNO2
                            % break the CNO sessions
                            % 2nd half of CNO3
                            if ~strcmpi(q5,badQ)
                                if RAND == true
                                    posAve = extract.cols(currentData,labels,'p','cell num',clusterNums{randCell},'session',sessions{5});
                                    blockLength = floor(size(posAve(:,1),1)/numBlocks);
                                    spikes = extract.cols(currentData,labels,'s','cell num',clusterNums{randCell},'session',sessions{5});
                                    [s5 t5 x5 y5] = sessionBlocks(2,2,blockLength,posAve,spikes);
                                else
                                    posAve = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'session',sessions{5});
                                    blockLength = floor(size(posAve(:,1),1)/numBlocks);
                                    spikes = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'session',sessions{5});
                                    [s5 t5 x5 y5] = sessionBlocks(2,2,blockLength,posAve,spikes);
                                end
                            end
                            
                            % 1st half of CNO4
                            if ~strcmpi(q6,badQ)
                                if RAND == true
                                    posAve = extract.cols(currentData,labels,'p','cell num',clusterNums{randCell},'session',sessions{6});
                                    blockLength = floor(size(posAve(:,1),1)/numBlocks);
                                    spikes = extract.cols(currentData,labels,'s','cell num',clusterNums{randCell},'session',sessions{6});
                                    [s6 t6 x6 y6] = sessionBlocks(1,2,blockLength,posAve,spikes);
                                else
                                    posAve = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'session',sessions{6});
                                    blockLength = floor(size(posAve(:,1),1)/numBlocks);
                                    spikes = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'session',sessions{6});
                                    [s6 t6 x6 y6] = sessionBlocks(1,2,blockLength,posAve,spikes);
                                end
                            end
                            
                            if ~strcmpi(q5,badQ) && ~strcmpi(q6,badQ)
                                offset = mean(diff(t5)) - min(t6) + t5(end);
                                sCNOday2 = [s5; s6+offset];
                                pCNOday2 = [t5 x5 y5; t6+offset x6 y6];
                                
                                map = analyses.map(pCNOday2,sCNOday2,'smooth',hippoGlobe.smoothing,'binWidth',hippoGlobe.binWidth,'limits',hippoGlobe.mapLimits);
                                map = map.z;
                                mapStore(:,:,iCluster,4,iGroup) = map;
                                
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
                                    rateStore(iCluster,iField,4,iGroup) = sortedRates(iField);
                                end
                                
                                if PLOT == true
                                    subplot(354)
                                    colorMapBRK(map);
                                    hold on
                                    for iField = 1:length(fieldStruct)
                                        rectangle('position',fieldStruct(iField).bbox,'linew',3)
                                    end
                                    subplot(359)
                                    bar(1:length(rates),rateStore(iCluster,1:length(rates),4,iGroup))
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    %% clean up data
    if iDataSet == 1
        % remove each group from big data store to make indexing easier
        eval('conAld = rateStore(:,:,:,1);');
        eval('expAld = rateStore(:,:,:,2);');
        
        rCheck = [];
        for i = 1:size(conAld,3)
            rCheck(i,:) = all(arrayfun(@isnan,conAld(:,:,i)'));
        end
        conAld = conAld(~all(rCheck),:,:);
        
        rCheck = [];
        for i = 1:size(expAld,3)
            rCheck(i,:) = all(arrayfun(@isnan,expAld(:,:,i)'));
        end
        expAld = expAld(~all(rCheck),:,:);
        
    elseif iDataSet == 2
        % remove each group from big data store to make indexing easier
        eval('conKA3 = rateStore(:,:,:,1);');
        
        rCheck = [];
        for i = 1:size(conKA3,3)
            rCheck(i,:) = all(arrayfun(@isnan,conKA3(:,:,i)'));
        end
        conKA3 = conKA3(~all(rCheck),:,:);
        
    elseif iDataSet == 3
        % remove each group from big data store to make indexing easier
        eval('conMUA = rateStore(:,:,:,1);');
        eval('expMUA = rateStore(:,:,:,2);');
        
        rCheck = [];
        for i = 1:size(conMUA,3)
            rCheck(i,:) = all(arrayfun(@isnan,conMUA(:,:,i)'));
        end
        conMUA = conMUA(~all(rCheck),:,:);
        
        rCheck = [];
        for i = 1:size(expMUA,3)
            rCheck(i,:) = all(arrayfun(@isnan,expMUA(:,:,i)'));
        end
        expMUA = expMUA(~all(rCheck),:,:);
    end
end

%% pool groups
conAll = vertcat(conAld,conKA3,conMUA);
expAll = vertcat(expAld,expMUA);

%% pool everything for random control
if RAND == true
    randAll = vertcat(conAll,expAll);
    expAll = randAll;
end

%% rate difference scatter plot - hM3
ratesBL1 = reshape(expAll(:,:,1),[],1);
ratesCNO1 = reshape(expAll(:,:,2),[],1);

diffDay1 = (ratesCNO1 - ratesBL1) ./ (ratesCNO1 + ratesBL1);

ratesBL2 = reshape(expAll(:,:,3),[],1);
ratesCNO2 = reshape(expAll(:,:,4),[],1);

diffDay2 = (ratesCNO2 - ratesBL2) ./ (ratesCNO2 + ratesBL2);

x = reshape(diffDay1,[],1);
y = reshape(diffDay2,[],1);
xnan = isnan(x);
ynan = isnan(y);
x = x(~xnan & ~ynan);
y = y(~xnan & ~ynan);

figure;
set(gcf,'color','w')
scatter(x,y,30,'k','filled')
hold on

r2 = calc.benFit(x,y,'degree',1,'showLine',1)
[r,p] = corr(x,y)
sum(~isnan(x))

set(gca,'box','off','fontsize',14,'fontweight','bold')
xlabel('Rate change BL-CNO1','fontsize',14,'fontweight','bold')
ylabel('Rate change BL-CNO2','fontsize',14,'fontweight','bold')
xlim([-1 1])
ylim([-1 1])
text(-0.95,0.95,sprintf('r = %.2f',r),'fontsize',14,'fontweight','bold')

%% grid field RD exp v con
compSesh = [1 2; 3 4; 1 3; 2 4];

expDiff = nan(size(expAll,1)*10,4);
for iComp = 1:length(compSesh)
    x = reshape(expAll(:,:,compSesh(iComp,1)),[],1);
    y = reshape(expAll(:,:,compSesh(iComp,2)),[],1);

    expDiff(:,iComp) = abs((x - y) ./ (x + y));
end

conDiff = nan(size(conAll,1)*10,4);
for iComp = 1:length(compSesh)

    x = reshape(conAll(:,:,compSesh(iComp,1)),[],1);
    y = reshape(conAll(:,:,compSesh(iComp,2)),[],1);

%     conDiff(:,iComp) = (y - x) ./ (y + x);
    conDiff(:,iComp) = abs((x - y) ./ (x + y));
end

nanmedian(expDiff)
nanmedian(conDiff)

%% bar graph + plot spread
figure;
set(gcf,'color','w')
hold on
xVals = [0.5 1.25 2.5 3.25];

h = bar(xVals(1:2:end),nanmedian(conDiff(:,1:2)));
set(h,'barwidth',0.25,'facecolor',[1 1 1],'edgecolor',[.8 .8 .8],'LineWidth',2.5)

h = bar(xVals(2:2:end),nanmedian(expDiff(:,1:2)));
set(h,'barwidth',0.25,'facecolor',[1 1 1],'edgecolor',[.5 .5 .5],'LineWidth',2.5)

plotSpread({conDiff(:,1) expDiff(:,1) conDiff(:,2) expDiff(:,2)},'xvalues',xVals);
o = flipud(findobj(gca,'type','line'));
set(o([1 3]),'markersize',8,'color',[0.8 0.8 0.8])
set(o([2 4]),'markersize',8,'color',[0.5 0.5 0.5])

set(gca,'xtickmode','manual','xtick',[0.87 2.87],'xticklabels',{'BL-CNO 1','BL-CNO 2'},'box','off','fontsize',14,'fontweight','bold')
ylabel('Grid field rate change','fontsize',14,'fontweight','bold')
xlim([-0.25 4])

%%
figure;
set(gcf,'color','w')
hold on
xVals = [0.5 1.25 2.5 3.25 4.5 5.25];

h = bar(xVals(1:2:end),nanmedian(conDiff(:,[3 1 2])));
set(h,'barwidth',0.25,'facecolor',[1 1 1],'edgecolor',[.8 .8 .8],'LineWidth',2.5)

h = bar(xVals(2:2:end),nanmedian(expDiff(:,[3 1 2])));
set(h,'barwidth',0.25,'facecolor',[1 1 1],'edgecolor',[.5 .5 .5],'LineWidth',2.5)

plotSpread({conDiff(:,3) expDiff(:,3) conDiff(:,1) expDiff(:,1) conDiff(:,2) expDiff(:,2)},'xvalues',xVals);
o = flipud(findobj(gca,'type','line'));
set(o([1 3 5]),'markersize',8,'color',[0.8 0.8 0.8])
set(o([2 4 6]),'markersize',8,'color',[0.5 0.5 0.5])

set(gca,'xtickmode','manual','xtick',[0.87 2.87 4.87],'xticklabels',{'BL-BL','BL-CNO 1','BL-CNO 2'},'box','off','fontsize',14,'fontweight','bold')
ylabel('Grid field rate change','fontsize',14,'fontweight','bold')
xlim([-0.25 6])
