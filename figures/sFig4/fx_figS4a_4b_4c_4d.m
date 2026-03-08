
global hippoGlobe
CYL 

% for S4a % S4c:
% measure = 'mean rate';

% for S4b % S4d:
measure = 'mean field size (cm2)';

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
thresh.spatInfo = 0.25;

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
toRemove = unique(extract.cols(raw,labels,'cell num','session',sessions{1},'mean rate','>=',thresh.FR));
raw = extract.rows(raw,labels,'remove','cell num',toRemove);

%% remove unstable cells
toRemove = unique(extract.cols(raw,labels,'cell num','session',sessions{1},'cc 1 vs 1','<',0.5));
raw = extract.rows(raw,labels,'remove','cell num',toRemove);

toRemove = unique(extract.cols(raw,labels,'cell num','session',sessions{1},'cc 1 vs 6','<',0.26));
raw = extract.rows(raw,labels,'remove','cell num',toRemove);

%% subset by group in case you want to plot something now
groupNames = {'con1x','conSal','exp1x'};
numGroups = length(groupNames);
con1x = extract.cols(raw,labels,':','group','con','injection','1x');
conSal = extract.cols(raw,labels,':','group','con','injection','sal');
exp1x = extract.cols(raw,labels,':','group','exp','injection','1x');

%% get values
% rate / size
compHalf = [1 2; 1 1; 2 2; 1 2];
compSesh = [1 1; 1 3; 1 3; 3 3];

% store data in 3 dimensions as: cell, session comparison, experimental group
dataStore = nan(1000,length(compSesh),numGroups);
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
                    % make sure quality isn't bad in either session
                    if ~strcmpi(q1,badQ) && ~strcmpi(q2,badQ)
                        % get values
                        rate1 = extract.cols(currentData,labels,measure,'cell num',clusterNums{iCluster},'exptHalf',expHalves{compHalf(iComp,1)},'session',sessions{compSesh(iComp,1)});
                        rate2 = extract.cols(currentData,labels,measure,'cell num',clusterNums{iCluster},'exptHalf',expHalves{compHalf(iComp,2)},'session',sessions{compSesh(iComp,2)});

                        % and store it
                        dataStore(iCluster,iComp,iGroup) = ((rate2 - rate1) / (rate2 + rate1));
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

%% BL-CNO only (bar graph + plot spread)
exp1xAbs = abs(exp1x);
conAllAbs = abs(conAll);

figure;
set(gcf,'color','w')
hold on
xVals = [0.5 1.25 2.5 3.25];

h = bar(xVals(1:2:end),nanmedian(conAllAbs(:,2:3)));
% set(h,'barwidth',0.3,'facecolor',[.8 .8 .8])
set(h,'barwidth',0.25,'facecolor',[1 1 1],'edgecolor',[.8 .8 .8],'LineWidth',2.5)

h = bar(xVals(2:2:end),nanmedian(exp1xAbs(:,2:3)));
% set(h,'barwidth',0.3,'facecolor',[.5 .5 .5])
set(h,'barwidth',0.25,'facecolor',[1 1 1],'edgecolor',[.5 .5 .5],'LineWidth',2.5)

plotSpread({conAllAbs(:,2) exp1xAbs(:,2) conAllAbs(:,3) exp1xAbs(:,3)},'xvalues',xVals);
o = flipud(findobj(gca,'type','line'));
set(o([1 3]),'markersize',8,'color',[0.8 0.8 0.8])
set(o([2 4]),'markersize',8,'color',[0.5 0.5 0.5])

set(gca,'xtickmode','manual','xtick',[0.87 2.87],'xticklabels',{'BL-CNO 1','BL-CNO 2'},'box','off','fontsize',14,'fontweight','bold')
% ylabel('Size change','fontsize',14,'fontweight','bold')
xlim([-0.25 4])

[p,~,stats] = ranksum(conAllAbs(:,2),exp1xAbs(:,2),'tail','left')
sum(~isnan(exp1xAbs(:,2)))
sum(~isnan(conAllAbs(:,2)))

[p,~,stats] = ranksum(conAllAbs(:,3),exp1xAbs(:,3),'tail','left')
sum(~isnan(exp1xAbs(:,3)))
sum(~isnan(conAllAbs(:,3)))

%% scatter plot of difference scores 
xVals = exp1x(:,2);
yVals = exp1x(:,3);
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

figure
set(gcf,'color','w')
plot(xVals,yVals,'k.','markersize',15)
hold on
axis([-1 1 -1 1])

[r,p] = corr(xVals,yVals)
r2 = calc.benFit(xVals,yVals,'degree',1,'showLine',1)

set(gca,'box','off','fontsize',14,'fontweight','bold')
xlabel('Size change BL-CNO 1','fontsize',14,'fontweight','bold')
ylabel('Size change BL-CNO 2','fontsize',14,'fontweight','bold')
text(-0.95,0.98,sprintf('r = %.2f',r),'fontsize',14,'fontweight','bold')

