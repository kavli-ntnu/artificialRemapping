
function raw = cleanUpQualityMEC(raw,labels,expHalves,numExpHalves,sessions,numSesh,badQ,offQ)

%% switch erroneous offQ to badQ
temp = {};
for iSession = 1:numSesh
    for iExpHalf = 1:numExpHalves
        currentData = extract.cols(raw,labels,':','session',sessions{iSession},'exptHalf',expHalves{iExpHalf});
        qToChange = extract.cols(currentData,labels,'cell num','quality',offQ,'mean rate','>=',0.1);
        if ~iscell(qToChange)
            qToChange = {qToChange};
        end
        for iCluster = 1:length(qToChange)
            currentData(strcmpi(currentData(:,strcmpi('cell num',labels)),qToChange{iCluster}) ...
                & strcmpi(currentData(:,strcmpi('exptHalf',labels)),expHalves{iExpHalf}), ...
                strcmpi('quality',labels)) = {badQ};
        end
        temp = [temp; currentData];
    end
end
raw = temp;

%% remove cells that are never good quality (Q3/Q4) in Expt Half 1
clusterNums = unique(extract.cols(raw,labels,'cell num'),'stable');
numClusters = length(clusterNums);
for iCluster = 1:numClusters
    quality1 = extract.cols(raw,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{1});
    quality = [quality1];
    if (sum(strcmpi(quality,offQ)) + sum(strcmpi(quality,badQ))) == length(quality)
        raw = extract.rows(raw,labels,'remove','cell num',clusterNums{iCluster});
    end
end

%% remove cells that are never good quality (Q3/Q4) in Expt Half 2
for iCluster = 1:numClusters
    quality2 = extract.cols(raw,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{2});
    quality = [quality2];
    if (sum(strcmpi(quality,offQ)) + sum(strcmpi(quality,badQ))) == length(quality)
        raw = extract.rows(raw,labels,'remove','cell num',clusterNums{iCluster});
    end
end

%% remove bad quality in either BL session
toRemove = extract.cols(raw,labels,'cell num','session',sessions{1},'quality',badQ,'exptHalf',expHalves{1});
if ~iscell(toRemove)
    toRemove = {toRemove};
end
raw = extract.rows(raw,labels,'remove','cell num',toRemove);
toRemove = extract.cols(raw,labels,'cell num','session',sessions{1},'quality',badQ,'exptHalf',expHalves{2});
if ~iscell(toRemove)
    toRemove = {toRemove};
end
raw = extract.rows(raw,labels,'remove','cell num',toRemove);
