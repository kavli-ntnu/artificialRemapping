
function fx_figS6b

global hippoGlobe
CYL

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
        
        doseNames = unique(extract.cols(raw,labels,'dose'),'sorted');
        
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
    
    %% subset by group in case you want to plot something now
    if iDataSet == 1
        groupNames = {'con1x','conSal','exp1x'};
        numGroups = length(groupNames);
        con1x = extract.cols(raw,labels,':','group','con','dose','cno');
        conSal = extract.cols(raw,labels,':','group','con','dose','veh');
        exp1x = extract.cols(raw,labels,':','group','exp','dose','1x');
        
    elseif iDataSet == 2
        groupNames = {'con1x','conSal','exp1x'};
        numGroups = length(groupNames);
        con1x = extract.cols(raw,labels,':','group','con','injection','1x');
        conSal = extract.cols(raw,labels,':','group','con','injection','sal');
        exp1x = extract.cols(raw,labels,':','group','exp','injection','1x');
        
    elseif iDataSet == 3
        groupNames = {'con1x','conSal','exp1x'};
        con1x = extract.cols(raw,labels,':','genotype','con','injection','CNO');
        conSal = extract.cols(raw,labels,':','genotype','exp','injection','sal');
        exp1x = extract.cols(raw,labels,':','genotype','exp','injection','CNO');
        
    elseif iDataSet == 4
        groupNames = {'con1x','conSal','exp1x'};
        con1x = extract.cols(raw,labels,':','genotype','con','injection','CNO');
        conSal = extract.cols(raw,labels,':','genotype','exp','injection','sal');
        exp1x = extract.cols(raw,labels,':','genotype','exp','injection','CNO');
    end
    
    %% get cc values
    % store data in 3 dimensions as: cell, session comparison, experimental group
    dataStore = nan(1000,3,numGroups);
    for iGroup = 1:numGroups
        eval(sprintf('currentData = %s;',groupNames{iGroup}))
        clusterNums = unique(extract.cols(currentData,labels,'cell num'),'stable');
        numClusters = length(clusterNums);
        for iCluster = 1:numClusters
            BLcheck = nan; CNOcheck = nan;
            if iDataSet == 1 || iDataSet == 2
                q1 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{1});
                q2 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{3});
                
                % must be a place cell in BL or CNO
                checkRateBL = extract.cols(currentData,labels,'peak rate','cell num',clusterNums{iCluster},'session',sessions{1}) > 0.5;
                checkFieldBL = extract.cols(currentData,labels,'number of fields','cell num',clusterNums{iCluster},'session',sessions{1}) > 0;
                checkInfoBL = extract.cols(currentData,labels,'spatial info','cell num',clusterNums{iCluster},'session',sessions{1}) > 0.67;     
                BLcheck = checkRateBL + checkFieldBL + checkInfoBL;
                
                checkRateCNO = extract.cols(currentData,labels,'peak rate','cell num',clusterNums{iCluster},'session',sessions{3}) > 0.5;
                checkFieldCNO = extract.cols(currentData,labels,'number of fields','cell num',clusterNums{iCluster},'session',sessions{3}) > 0;
                checkInfoCNO = extract.cols(currentData,labels,'spatial info','cell num',clusterNums{iCluster},'session',sessions{3}) > 0.67;     
                CNOcheck = checkRateCNO + checkFieldCNO + checkInfoCNO;
                
                % make sure quality isn't bad in either session && ~strcmpi(q1,offQ) && ~strcmpi(q2,offQ)
                if ~strcmpi(q1,badQ) && ~strcmpi(q2,badQ)
                    if BLcheck == 3 || CNOcheck == 3
                        % get cc value
                        map1 = extract.cols(currentData,labels,'rate map','cell num',clusterNums{iCluster},'session',sessions{1});
                        map2 = extract.cols(currentData,labels,'rate map','cell num',clusterNums{iCluster},'session',sessions{3});
                        cc = analyses.spatialCrossCorrelation(map1,map2);
                        
                        % and store it
                        dataStore(iCluster,1,iGroup) = str2num(clusterNums{iCluster});
                        dataStore(iCluster,2,iGroup) = extract.cols(currentData,labels,'CC 1 vs 1','cell num',clusterNums{iCluster},'session',sessions{1});
                        dataStore(iCluster,3,iGroup) = cc;
                    end
                end
                
            elseif iDataSet == 3 || iDataSet == 4
                q1 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{1});
                q2 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'session',sessions{2});
                
                % make sure quality isn't bad in either session && ~strcmpi(q1,offQ) && ~strcmpi(q2,offQ)
                if ~strcmpi(q1,badQ) && ~strcmpi(q2,badQ)
                    numBlocks = 2;
                    
                    map1 = extract.cols(currentData,labels,'rate map','cell num',clusterNums{iCluster},'session',sessions{1});
                    %
                    if ~strcmpi(q2,badQ)
                        posAve = extract.cols(currentData,labels,'p','cell num',clusterNums{iCluster},'session',sessions{2});
                        blockLength = floor(size(posAve(:,1),1)/numBlocks);
                        spikes = extract.cols(currentData,labels,'s','cell num',clusterNums{iCluster},'session',sessions{2});
                        [s2 t2 x2 y2] = sessionBlocks(2,numBlocks,blockLength,posAve,spikes);
                        map2 = analyses.map([t2 x2 y2],s2,'smooth',hippoGlobe.smoothing,'binWidth',hippoGlobe.binWidth,'limits',hippoGlobe.mapLimits);
                        
                        [~,CNOfield] = analyses.placefield(map2.z);

                        % must be a place cell in BL or CNO
                        checkRateBL = extract.cols(currentData,labels,'peak rate','cell num',clusterNums{iCluster},'session',sessions{1}) > 0.5;
                        checkFieldBL = extract.cols(currentData,labels,'number of fields','cell num',clusterNums{iCluster},'session',sessions{1}) > 0;
                        checkInfoBL = extract.cols(currentData,labels,'spatial info','cell num',clusterNums{iCluster},'session',sessions{1}) > 0.67;
                        BLcheck = checkRateBL + checkFieldBL + checkInfoBL;
                        
                        checkFieldCNO = length(CNOfield) > 0;
%                         display(iGroup)
%                         display(iCluster)
                        if checkFieldCNO == 1
                            checkRateCNO = nanmax([CNOfield.peak]) > 0.5;
                        else
                            checkRateCNO = 0;
                        end
                        checkInfoCNO = extract.cols(currentData,labels,'spatial info','cell num',clusterNums{iCluster},'session',sessions{2}) > 0.67;
                        CNOcheck = checkRateCNO + checkFieldCNO + checkInfoCNO;
                        
                        if BLcheck == 3 || CNOcheck == 3
                            dataStore(iCluster,1,iGroup) = str2num(clusterNums{iCluster});
                            dataStore(iCluster,2,iGroup) = extract.cols(currentData,labels,'CC 1 vs 1','cell num',clusterNums{iCluster},'session',sessions{1});
                        
                            % get cc values
                            cc = analyses.spatialCrossCorrelation(map1,map2.z);
                            dataStore(iCluster,3,iGroup) = cc;
                        end                   
                    end
                end
            end
        end
    end


    %% clean up data
    if iDataSet == 1
        % remove each group from big data store to make indexing easier
        eval('con1xAld = dataStore(:,:,1);');
        eval('conSalAld = dataStore(:,:,2);');
        eval('exp1xAld = dataStore(:,:,3);');
        
        % remove rows of all nans that were used to pad the array
        con1xAld(all(arrayfun(@isnan,con1xAld)'),:) = [];
        conSalAld(all(arrayfun(@isnan,conSalAld)'),:) = [];
        exp1xAld(all(arrayfun(@isnan,exp1xAld)'),:) = [];
        
    elseif iDataSet == 2
        % remove each group from big data store to make indexing easier
        eval('con1xJas = dataStore(:,:,1);');
        eval('conSalJas = dataStore(:,:,2);');
        eval('exp1xJas = dataStore(:,:,3);');
        
        % remove rows of all nans that were used to pad the array
        con1xJas(all(arrayfun(@isnan,con1xJas)'),:) = [];
        conSalJas(all(arrayfun(@isnan,conSalJas)'),:) = [];
        exp1xJas(all(arrayfun(@isnan,exp1xJas)'),:) = [];
        
    elseif iDataSet == 3
        % remove each group from big data store to make indexing easier
        eval('con1xCont = dataStore(:,:,1);');
        eval('conSalCont = dataStore(:,:,2);');
        eval('exp1xCont = dataStore(:,:,3);');
        
        % remove rows of all nans that were used to pad the array
        con1xCont(all(arrayfun(@isnan,con1xCont)'),:) = [];
        conSalCont(all(arrayfun(@isnan,conSalCont)'),:) = [];
        exp1xCont(all(arrayfun(@isnan,exp1xCont)'),:) = [];
        
    elseif iDataSet == 4
        % remove each group from big data store to make indexing easier
        eval('con1xBreak = dataStore(:,:,1);');
        eval('conSalBreak = dataStore(:,:,2);');
        eval('exp1xBreak = dataStore(:,:,3);');
        
        % remove rows of all nans that were used to pad the array
        con1xBreak(all(arrayfun(@isnan,con1xBreak)'),:) = [];
        conSalBreak(all(arrayfun(@isnan,conSalBreak)'),:) = [];
        exp1xBreak(all(arrayfun(@isnan,exp1xBreak)'),:) = [];
    end
end

% pool groups
exp1x = vertcat(exp1xJas,exp1xAld,exp1xCont,exp1xBreak);
nExp1x = num2str(max(sum(~isnan(exp1x))))

conAll = vertcat(con1xJas,conSalJas,con1xAld,conSalAld,con1xCont,conSalCont,con1xBreak,conSalBreak);
nConAll = num2str(max(sum(~isnan(conAll))))

nanmedian(exp1x)
nanmedian(conAll)

sum(~isnan(exp1x))
sum(~isnan(conAll))

hM3 = exp1x(:,3);

%% compare to degree of remapping in simulation
% grid subfield rate change, distribution median = 0.6
load fig4b

spatialCorrAll = nan(200*numSims,1);
for i = 1:200*numSims
    if iBL(i) || iCNO(i)
        
        mapBL = general.smooth(pBL(:,:,i),6.67);
        mapCNO = general.smooth(pCNO(:,:,i),6.67);
        spatialCorrAll(i) = analyses.spatialCrossCorrelation(mapBL,mapCNO);
    end
end

nanmedian(spatialCorrAll)
sum(~isnan(spatialCorrAll))

sc6 = spatialCorrAll;

%% load data
% grid subfield rate change, distribution median = 0.3
load figS6b_gfc3

spatialCorrAll = nan(200*numSims,1);
for i = 1:200*numSims
    if iBL(i) || iCNO(i)
        mapBL = general.smooth(pBL(:,:,i),6.67);
        mapCNO = general.smooth(pCNO(:,:,i),6.67);
        spatialCorrAll(i) = analyses.spatialCrossCorrelation(mapBL,mapCNO);
    end
end

nanmedian(spatialCorrAll)
sum(~isnan(spatialCorrAll))

sc3 = spatialCorrAll;

%% load data
% grid subfield rate change, distribution median = 0.9
load figS6b_gfc9

spatialCorrAll = nan(200*numSims,1);
for i = 1:200*numSims
    if iBL(i) || iCNO(i)
        mapBL = general.smooth(pBL(:,:,i),6.67);
        mapCNO = general.smooth(pCNO(:,:,i),6.67);
        spatialCorrAll(i) = analyses.spatialCrossCorrelation(mapBL,mapCNO);
    end
end

nanmedian(spatialCorrAll)

sum(~isnan(spatialCorrAll))

sc9 = spatialCorrAll;

%%
% data plotted in figS6b