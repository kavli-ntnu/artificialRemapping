
function raw = cleanUpQualityTilRom(raw,labels,sessions,numSesh,badQ,offQ)

%% switch erroneous offQ to badQ
temp = {};
for iSession = 1:numSesh
    currentData = extract.cols(raw,labels,':','session',sessions{iSession});
    qToChange = extract.cols(currentData,labels,'cell num','quality',offQ,'mean rate','>=',0.15);
    if ~iscell(qToChange)
        qToChange = {qToChange};
    end
    for iCluster = 1:length(qToChange)
        currentData(strcmpi(currentData(:,strcmpi('cell num',labels)),qToChange{iCluster}), ...
            strcmpi('quality',labels)) = {badQ};
    end
    temp = [temp; currentData];
end
raw = temp;

%% remove cells that are never good quality
clusterNums = unique(extract.cols(raw,labels,'cell num'),'stable');
numClusters = length(clusterNums);
for iCluster = 1:numClusters
    quality = extract.cols(raw,labels,'quality','cell num',clusterNums{iCluster});
    if (sum(strcmpi(quality,offQ)) + sum(strcmpi(quality,badQ))) == length(quality)
        raw = extract.rows(raw,labels,'remove','cell num',clusterNums{iCluster});
    end
end

%% remove bad quality in either BL session
toRemove = extract.cols(raw,labels,'cell num','session',sessions{1},'quality',badQ);
if ~iscell(toRemove)
    toRemove = {toRemove};
end
raw = extract.rows(raw,labels,'remove','cell num',toRemove);

toRemove = extract.cols(raw,labels,'cell num','session',sessions{2},'quality',badQ);
if ~iscell(toRemove)
    toRemove = {toRemove};
end
raw = extract.rows(raw,labels,'remove','cell num',toRemove);
