
function fx_fig2b_2c

global hippoGlobe
CYL

%% load data
load masterMatJasmineU2_new
raw = dataOutput;

%% switch numbers to strings so you can use functions like 'unique'
for iRow = 1:size(raw,1)
    raw{iRow,strcmpi('cell num',labels)} = num2str(raw{iRow,strcmpi('cell num',labels)});
    raw{iRow,strcmpi('quality',labels)} = num2str(raw{iRow,strcmpi('quality',labels)});
    raw{iRow,strcmpi('exptHalf',labels)} = num2str(raw{iRow,strcmpi('exptHalf',labels)});
end

%% functional thresholds
thresh.FR = 7;
thresh.spatInfo = 1.5;

%% store names of groups and sessions
injectionNames = unique(extract.cols(raw,labels,'injection'),'sorted');
numInjectionGroups = numel(injectionNames);
sessions = unique(extract.cols(raw,labels,'session'),'stable');
numSesh = numel(sessions);
expHalves = unique(extract.cols(raw,labels,'exptHalf'),'stable');
numExpHalves = numel(expHalves);

%% fix quality
badQ = '3';
offQ = '4';
% see documentation for this function to see how it cleans the data
raw = cleanUpQualityCML(raw,labels,expHalves,numExpHalves,sessions,numSesh,badQ,offQ);

%% only look at excitatory cells
toRemove = unique(extract.cols(raw,labels,'cell num','session',sessions{1},'mean rate','>',thresh.FR));
% now we remove all those inhibitory cells we just found
raw = extract.rows(raw,labels,'remove','cell num',toRemove);

%% remove unstable cells
toRemove = unique(extract.cols(raw,labels,'cell num','session',sessions{1},'cc 1 vs 1','<',0.5));
raw = extract.rows(raw,labels,'remove','cell num',toRemove);

toRemove = unique(extract.cols(raw,labels,'cell num','session',sessions{1},'cc 1 vs 6','<',0.26));
raw = extract.rows(raw,labels,'remove','cell num',toRemove);

%% subset by group in case you want to plot something now
% note we use the colon to get all columns instead of just one
groupNames = {'con1x','conSal','exp1x'};
numGroups = length(groupNames);
con1x = extract.cols(raw,labels,':','group','con','injection','1x');
conSal = extract.cols(raw,labels,':','group','con','injection','sal');
exp1x = extract.cols(raw,labels,':','group','exp','injection','1x');

%% get cc values
compHalf = [1 2; 1 1; 2 2; 1 2];
compSesh = [1 1; 1 5; 1 5; 5 3];

% store data in 3 dimensions as: cell, session comparison, experimental group
dataStore = nan(1000,length(compSesh)+1,numGroups);
for iGroup = 1:numGroups
    eval(sprintf('currentData = %s;',groupNames{iGroup}))
    clusterNums = unique(extract.cols(currentData,labels,'cell num'),'stable');
    numClusters = length(clusterNums);
    for iCluster = 1:numClusters
        mapBL1 = extract.cols(currentData,labels,'rate map','cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{1});
        
        [~,fieldsBL1] = analyses.placefield(mapBL1,'threshold',0.4,'binWidth',hippoGlobe.binWidth,'minBins',10,'minPeak',0.5);
        qBL1 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{1});
        
        mapBL2 = extract.cols(currentData,labels,'rate map','cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{1});
        [~,fieldsBL2] = analyses.placefield(mapBL2,'threshold',0.4,'binWidth',hippoGlobe.binWidth,'minBins',10,'minPeak',0.5);
        qBL2 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{2},'session',sessions{1});
        
        if ~isempty(fieldsBL1) && ~strcmpi(qBL1,badQ) && ~isempty(fieldsBL2) && ~strcmpi(qBL2,badQ)
            for iComp = 1:length(compSesh)
                try     % use a try block to keep going if there are missing values
                    % get quality for each session to compare
                    q1 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{compHalf(iComp,1)},'session',sessions{compSesh(iComp,1)});
                    q2 = extract.cols(currentData,labels,'quality','cell num',clusterNums{iCluster},'exptHalf',expHalves{compHalf(iComp,2)},'session',sessions{compSesh(iComp,2)});
                    
                    % make sure quality isn't bad in either session && ~strcmpi(q1,offQ) && ~strcmpi(q2,offQ)
                    if ~strcmpi(q1,badQ) && ~strcmpi(q2,badQ)
                        % get cc value
                        cc = extract.cols(currentData,labels,sprintf('cc %d vs %d',compSesh(iComp,1)+(5*(compHalf(iComp,1)-1)),compSesh(iComp,2)+(5*(compHalf(iComp,2)-1))),'cell num',clusterNums{iCluster},'exptHalf',expHalves{1},'session',sessions{1}); % just need any session here
                        % and store it
                        dataStore(iCluster,1,iGroup) = str2num(clusterNums{iCluster});
                        dataStore(iCluster,iComp+1,iGroup) = cc;
                    end
                end
            end
        end
    end
end

%% clean up data
% remove each group from big data store to make indexing easier
eval('con1x = dataStore(:,:,1);');
eval('conSal = dataStore(:,:,2);');
eval('exp1x = dataStore(:,:,3);');

% remove rows of all nans that were used to pad the array
con1x(all(arrayfun(@isnan,con1x)'),:) = [];
conSal(all(arrayfun(@isnan,conSal)'),:) = [];
exp1x(all(arrayfun(@isnan,exp1x)'),:) = [];

% count the max number of cells for any session
nCon1x = num2str(max(sum(~isnan(con1x))));
nConSal = num2str(max(sum(~isnan(conSal))));
nExp1x = num2str(max(sum(~isnan(exp1x))));

%% pool control groups
conAll = vertcat(con1x,conSal);
nConAll = num2str(max(sum(~isnan(conAll))));

nanmedian(exp1x)
sum(~isnan(exp1x))

%% BLxCNO1 & BLxCNO2
figure;
[f,x,flo,fup] = ecdf(conAll(:,3));
[h1,h2,h3] = shadedplot(x',flo',fup',[219/255 61/255 49/255]);
set(h1,'FaceAlpha',0.2,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(conAll(:,3));
set(C,'color',[219/255 61/255 49/255 0.4],'linewidth',3)

hold on
[f,x,flo,fup] = ecdf(conAll(:,4));
[h1,h2,h3] = shadedplot(x',flo',fup',[219/255 61/255 49/255]);
set(h1,'FaceAlpha',0.3,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(conAll(:,4));
set(C,'color',[219/255 61/255 49/255],'linewidth',3)

[f,x,flo,fup] = ecdf(exp1x(:,3));
[h1,h2,h3] = shadedplot(x',flo',fup',[10/255 113/255 178/255]);
set(h1,'FaceAlpha',0.2,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(exp1x(:,3));
set(C,'color',[10/255 113/255 178/255 0.4],'linewidth',3)

hold on
[f,x,flo,fup] = ecdf(exp1x(:,4));
[h1,h2,h3] = shadedplot(x',flo',fup',[10/255 113/255 178/255]);
set(h1,'FaceAlpha',0.3,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(exp1x(:,4));
set(C,'color',[10/255 113/255 178/255],'linewidth',3)

xlabel('Spatial correlation')
ylabel('Cumulative proprotion')
title ''
grid off
set(gca,'ytick',0:0.25:1)

%% CNO1xCNO2
figure;
[f,x,flo,fup] = ecdf(conAll(:,5));
[h1,h2,h3] = shadedplot(x',flo',fup',[132/255 32/255 34/255]);
set(h1,'FaceAlpha',0.3,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(conAll(:,5));
set(C,'color',[132/255 32/255 34/255],'linewidth',3)

hold on
[f,x,flo,fup] = ecdf(exp1x(:,5));
[h1,h2,h3] = shadedplot(x',flo',fup',[25/255 72/255 101/255]);
set(h1,'FaceAlpha',0.3,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(exp1x(:,5));
set(C,'color',[25/255 72/255 101/255],'linewidth',3)

xlabel('Spatial correlation')
ylabel('Cumulative proprotion')
title ''
grid off
set(gca,'ytick',0:0.2:1)

%% BL1xBL2
figure;
[f,x,flo,fup] = ecdf(conAll(:,2));
[h1,h2,h3] = shadedplot(x',flo',fup',[252/255 186/255 131/255]);
set(h1,'FaceAlpha',0.3,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(conAll(:,2));
set(C,'color',[252/255 186/255 131/255],'linewidth',3)

hold on
[f,x,flo,fup] = ecdf(exp1x(:,2));
[h1,h2,h3] = shadedplot(x',flo',fup',[166/255 189/255 220/255]);
set(h1,'FaceAlpha',0.3,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(exp1x(:,2));
set(C,'color',[166/255 189/255 220/255],'linewidth',3)

xlabel('Spatial correlation')
ylabel('Cumulative proprotion')
title ''
grid off
set(gca,'ytick',0:0.25:1)

%%
sum(~isnan(exp1x))
sum(~isnan(conAll))

[~,p,ks] = kstest2(exp1x(:,3),conAll(:,3))
[~,p,ks] = kstest2(exp1x(:,4),conAll(:,4))
[~,p,ks] = kstest2(exp1x(:,5),conAll(:,5))
