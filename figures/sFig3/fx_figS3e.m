function fx_figS3e

global hippoGlobe

%% check inputs and load data
cellType = 'excitatory';
measure = 'rate map';
pkThresh = 1;

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
    
    mapStore = cell(1);
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
        
        % fix quality
        raw = cleanUpQualityMEC(raw,labels,expHalves,numExpHalves,sessions,numSesh,badQ,offQ);
    elseif iDataSet == 2
        % fix quality
        raw = extract.cleanUpQuality(raw,labels,sessions,numSesh,badQ,offQ);
    elseif iDataSet == 3
        % fix quality
        raw = extract.cleanUpQuality(raw,labels,sessions,numSesh,badQ,offQ);
    end
    
    %% functional thresholds
    fieldThresh = 0.35;
    fieldPeak = 0.1;
    minPeak = 1;
    
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
        
    elseif iDataSet == 2
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
            
        elseif iDataset == 3
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
    end
    
    %% separate functional groups
    if iDataSet == 1
        groupNames = unique(extract.cols(raw,labels,'group'),'sorted');
        numGroups = numel(groupNames);
        Con = raw(strcmpi(raw(:,strcmpi('group',labels)),'con'),:);
        Exp = raw(strcmpi(raw(:,strcmpi('group',labels)),'exp'),:);
        
    elseif iDataSet == 2
        groupNames = unique(extract.cols(raw,labels,'group'),'sorted');
        numGroups = numel(groupNames);
        Con = raw(strcmpi(raw(:,strcmpi('group',labels)),'con'),:);
        
    elseif iDataSet == 3
        groupNames = {'Con','Exp'};
        numGroups = numel(groupNames);
        Con = extract.cols(raw,labels,':','genotype','exp','injection','sal');
        Exp = extract.cols(raw,labels,':','genotype','exp','injection','cno');
    end
    
    %% initialize
    compSesh = [1 2; 3 4; 1 3; 2 4];
    
    if iDataSet == 1
        dataStore = nan(1000,length(compSesh),numGroups);
        
        for iGroup = 1:numGroups
            eval(sprintf('currentData = %s;',groupNames{iGroup}))
            clusterNums = unique(extract.cols(currentData,labels,'cell num'),'stable');
            numClusters = length(clusterNums);
            
            for iCluster = 1:numClusters
                mapTemp = cell(1,4);
                
                q1 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{1});
                q2 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{3});
                q3 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{4});
                q4 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{1});
                q5 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{3});
                q6 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{4});
                
                % rates for baselines
                if ~strcmpi(q1,badQ)
                    mapTemp{1,1} = extract.cols(currentData,labels,measure,'cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{1});
                end
                if ~strcmpi(q4,badQ)
                    mapTemp{1,3} = extract.cols(currentData,labels,measure,'cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{1});
                end
                
                % C1
                if ~strcmpi(q2,badQ)
                    p2 = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{3});
                    s2 = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{3});
                end
                
                % D1
                if ~strcmpi(q3,badQ)
                    p3 = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{4});
                    s3 = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{4});
                end
                
                % C2
                if ~strcmpi(q5,badQ)
                    p5 = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{3});
                    s5 = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{3});
                end
                
                % D2
                if ~strcmpi(q6,badQ)
                    p6 = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{4});
                    s6 = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{4});
                end
                
                % merge sessions
                if ~strcmpi(q2,badQ) && ~strcmpi(q3,badQ)
                    offset = mean(diff(p2(:,1))) - min(p3(:,1)) + p2(end,1);
                    sCNOday1 = [s2; s3+offset];
                    pCNOday1 = [p2; p3+offset];
                    map = analyses.map(pCNOday1,sCNOday1,'smooth',hippoGlobe.smoothing,'binWidth',hippoGlobe.binWidth,'limits',hippoGlobe.mapLimits);
                    mapTemp{1,2} = map.z;
                end
                if ~strcmpi(q5,badQ) && ~strcmpi(q6,badQ)
                    offset = mean(diff(p5(:,1))) - min(p6(:,1)) + p5(end,1);
                    sCNOday2 = [s5; s6+offset];
                    pCNOday2 = [p5; p6+offset];
                    map = analyses.map(pCNOday2,sCNOday2,'smooth',hippoGlobe.smoothing,'binWidth',hippoGlobe.binWidth,'limits',hippoGlobe.mapLimits);
                    mapTemp{1,4} = map.z;
                end
                
                sc = nan(1,4);
                for iComp = 1:size(compSesh,1)
                    if nanmax(nanmax(mapTemp{1,compSesh(iComp,1)})) > pkThresh & nanmax(nanmax(mapTemp{1,compSesh(iComp,2)})) > pkThresh
                        sc(1,iComp) = analyses.spatialCrossCorrelation(mapTemp{1,compSesh(iComp,1)},mapTemp{1,compSesh(iComp,2)});
                    end
                end
                
                dataStore(iCluster,1:size(compSesh,1),iGroup) = sc;
                mapStore{iCluster,:,iGroup} = mapTemp;
            end
        end
        
    elseif iDataSet == 2
        
        dataStore = nan(1000,length(compSesh),numGroups);
        
        for iGroup = 1:numGroups
            currentData = Con;
            clusterNums = unique(extract.cols(currentData,labels,'cell num'),'stable');
            numClusters = length(clusterNums);
            
            numBlocks = 4;
            
            for iCluster = 1:numClusters
                
                mapTemp = cell(1,4);
                
                % quality for each piece
                q1 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{1});
                q2 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{2});
                q3 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{3});
                q4 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{4});
                
                % rates for baselines
                if ~strcmpi(q1,badQ)
                    mapTemp{1,1} = extract.cols(currentData,labels,measure,'cell num',clusterNums{iCluster},'session',sessions{1});
                end
                if ~strcmpi(q3,badQ)
                    mapTemp{1,3} = extract.cols(currentData,labels,measure,'cell num',clusterNums{iCluster},'session',sessions{3});
                end
                
                % break the CNO sessions
                if ~strcmpi(q2,badQ)
                    posAve = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'session',sessions{2});
                    blockLength = floor(size(posAve(:,1),1)/numBlocks);
                    spikes = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'session',sessions{2});
                    [s2 t2 x2 y2] = sessionBlocks(2,numBlocks,blockLength,posAve,spikes);
                    [s3 t3 x3 y3] = sessionBlocks(3,numBlocks,blockLength,posAve,spikes);
                end
                
                if ~strcmpi(q4,badQ)
                    posAve = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'session',sessions{4});
                    blockLength = floor(size(posAve(:,1),1)/numBlocks);
                    spikes = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'session',sessions{4});
                    [s5 t5 x5 y5] = sessionBlocks(2,numBlocks,blockLength,posAve,spikes);
                    [s6 t6 x6 y6] = sessionBlocks(3,numBlocks,blockLength,posAve,spikes);
                end
                
                % merge blocks
                if ~strcmpi(q2,badQ)
                    offset = mean(diff(t2)) - min(t3) + t2(end);
                    sCNOday1 = [s2; s3+offset];
                    pCNOday1 = [t2 x2 y2; t3+offset x3 y3];
                    map = analyses.map(pCNOday1,sCNOday1,'smooth',hippoGlobe.smoothing,'binWidth',hippoGlobe.binWidth,'limits',hippoGlobe.mapLimits);
                    mapTemp{1,2} = map.z;
                end
                if ~strcmpi(q4,badQ)
                    offset = mean(diff(t5)) - min(t6) + t5(end);
                    sCNOday2 = [s5; s6+offset];
                    pCNOday2 = [t5 x5 y5; t6+offset x6 y6];
                    map = analyses.map(pCNOday2,sCNOday2,'smooth',hippoGlobe.smoothing,'binWidth',hippoGlobe.binWidth,'limits',hippoGlobe.mapLimits);
                    mapTemp{1,4} = map.z;
                end
                
                sc = nan(1,4);
                for iComp = 1:size(compSesh,1)
                    if nanmax(nanmax(mapTemp{1,compSesh(iComp,1)})) > pkThresh & nanmax(nanmax(mapTemp{1,compSesh(iComp,2)})) > pkThresh
                        sc(1,iComp) = analyses.spatialCrossCorrelation(mapTemp{1,compSesh(iComp,1)},mapTemp{1,compSesh(iComp,2)});
                    end
                end
                
                dataStore(iCluster,1:size(compSesh,1),iGroup) = sc;
                mapStore{iCluster,:,iGroup} = mapTemp;
            end
        end
        
        
        
    elseif iDataSet == 3
        
        dataStore = nan(1000,length(compSesh),numGroups);
        
        for iGroup = 1:numGroups
            eval(sprintf('currentData = %s;',groupNames{iGroup}))
            clusterNums = unique(extract.cols(currentData,labels,'cell num'),'stable');
            numClusters = length(clusterNums);
            
            numBlocks = 2;
            
            for iCluster = 1:numClusters
                
                mapTemp = cell(1,4);
                
                % quality for each piece
                q1 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{1});
                q2 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{2});
                q3 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{3});
                q4 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{4});
                q5 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{5});
                q6 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{6});
                
                % rates for baselines
                if ~strcmpi(q1,badQ)
                    mapTemp{1,1} = extract.cols(currentData,labels,measure,'cell num',clusterNums{iCluster},'session',sessions{1});
                end
                if ~strcmpi(q4,badQ)
                    mapTemp{1,3} = extract.cols(currentData,labels,measure,'cell num',clusterNums{iCluster},'session',sessions{4});
                end
                
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
                
                % 2nd half of CNO3
                if ~strcmpi(q5,badQ)
                    posAve = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'session',sessions{5});
                    blockLength = floor(size(posAve(:,1),1)/numBlocks);
                    spikes = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'session',sessions{5});
                    [s5 t5 x5 y5] = sessionBlocks(2,2,blockLength,posAve,spikes);
                end
                
                % 1st half of CNO4
                if ~strcmpi(q6,badQ)
                    posAve = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'session',sessions{6});
                    blockLength = floor(size(posAve(:,1),1)/numBlocks);
                    spikes = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'session',sessions{6});
                    [s6 t6 x6 y6] = sessionBlocks(1,2,blockLength,posAve,spikes);
                end
                
                % merge blocks
                if ~strcmpi(q2,badQ) && ~strcmpi(q3,badQ)
                    offset = mean(diff(t2)) - min(t3) + t2(end);
                    sCNOday1 = [s2; s3+offset];
                    pCNOday1 = [t2 x2 y2; t3+offset x3 y3];
                    map = analyses.map(pCNOday1,sCNOday1,'smooth',hippoGlobe.smoothing,'binWidth',hippoGlobe.binWidth,'limits',hippoGlobe.mapLimits);
                    mapTemp{1,2} = map.z;
                end
                if ~strcmpi(q5,badQ) && ~strcmpi(q6,badQ)
                    offset = mean(diff(t5)) - min(t6) + t5(end);
                    sCNOday2 = [s5; s6+offset];
                    pCNOday2 = [t5 x5 y5; t6+offset x6 y6];
                    map = analyses.map(pCNOday2,sCNOday2,'smooth',hippoGlobe.smoothing,'binWidth',hippoGlobe.binWidth,'limits',hippoGlobe.mapLimits);
                    mapTemp{1,4} = map.z;
                end
                
                sc = nan(1,4);
                for iComp = 1:size(compSesh,1)
                    if nanmax(nanmax(mapTemp{1,compSesh(iComp,1)})) > pkThresh & nanmax(nanmax(mapTemp{1,compSesh(iComp,2)})) > pkThresh
                        sc(1,iComp) = analyses.spatialCrossCorrelation(mapTemp{1,compSesh(iComp,1)},mapTemp{1,compSesh(iComp,2)});
                    end
                end
                
                dataStore(iCluster,1:size(compSesh,1),iGroup) = sc;
                mapStore{iCluster,:,iGroup} = mapTemp;
            end
        end
    end
    
    
    %% clean up data
    if iDataSet == 1
        % remove each group from big data store to make indexing easier
        eval('conAld = dataStore(:,:,1);');
        eval('expAld = dataStore(:,:,2);');
        
        conAld(all(arrayfun(@isnan,conAld)'),:) = [];
        expAld(all(arrayfun(@isnan,expAld)'),:) = [];
        
        conAldMaps = mapStore(:,:,1);
        expAldMaps = mapStore(:,:,2);
        
    elseif iDataSet == 2
        % remove each group from big data store to make indexing easier
        eval('conKA3 = dataStore(:,:,1);');
        
        conKA3(all(arrayfun(@isnan,conKA3)'),:) = [];
        
        conKA3Maps = mapStore(:,:,1);
        
    elseif iDataSet == 3
        % remove each group from big data store to make indexing easier
        eval('conMUA = dataStore(:,:,1);');
        eval('expMUA = dataStore(:,:,2);');
        
        conMUA(all(arrayfun(@isnan,conMUA)'),:) = [];
        expMUA(all(arrayfun(@isnan,expMUA)'),:) = [];
        
        conMUAMaps = mapStore(:,:,1);
        expMUAMaps = mapStore(:,:,2);
    end
end

conAll = vertcat(conAld,conKA3,conMUA);
expAll = vertcat(expAld,expMUA);
nanmedian(conAll)
nanmedian(expAll)

%%
conAllMaps = vertcat(conAldMaps,conKA3Maps,conMUAMaps);
conAllMaps = minions.removeEmpties(conAllMaps,'rows','all');

expAllMaps = vertcat(expAldMaps,expMUAMaps);
expAllMaps = minions.removeEmpties(expAllMaps,'rows','all');


%% bar graph + plot spread
figure;
set(gcf,'color','w')
hold on
xVals = [0.5 1.25 2.5 3.25];

plotSpread({conAll(:,1) expAll(:,1) conAll(:,2) expAll(:,2)},'xvalues',xVals);
o = flipud(findobj(gca,'type','line'));
set(o([1 3]),'markersize',8,'color',[0.8 0.8 0.8])
set(o([2 4]),'markersize',8,'color',[0.5 0.5 0.5])

set(gca,'xtickmode','manual','xtick',[0.87 2.87],'xticklabels',{'BL-CNO 1','BL-CNO 2'},'box','off','fontsize',14,'fontweight','bold')
ylabel('Spatial correlation','fontsize',14,'fontweight','bold')
xlim([-0.25 4])

%%
[p,~,stats] = ranksum(conAll(:,1),expAll(:,1))
sum(~isnan(conAll(:,1)))
sum(~isnan(expAll(:,1)))

[p,~,stats] = ranksum(conAll(:,2),expAll(:,2))
sum(~isnan(conAll(:,2)))
sum(~isnan(expAll(:,2)))
