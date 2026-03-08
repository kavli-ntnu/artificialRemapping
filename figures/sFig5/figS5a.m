startup
CYL
hippoGlobe.arena
load masterMat_HP2

%%
currentData = dataOutput;
sessions = unique(extract.cols(currentData,labels,'session'),'stable');
for iRow = 1:size(currentData,1)
    currentData{iRow,strcmpi('exptHalf',labels)} = num2str(currentData{iRow,strcmpi('exptHalf',labels)});
end
expHalves = unique(extract.cols(currentData,labels,'exptHalf'),'sorted');
clusterNums = unique(extract.cols(currentData,labels,'cell num'),'stable');

%%
selectedCells = [283 405 96 354 314 377 481 1139];

for i = 1:size(selectedCells,2)
    figure;
    cellNum = mat2str(selectedCells(i));
    
    posAve = extract.cols(currentData,labels,'p','cell num',cellNum,'session',sessions{1},'exptHalf',expHalves{1});
    spikes = extract.cols(currentData,labels,'s','cell num',cellNum,'session',sessions{1},'exptHalf',expHalves{1});
    spos = extract.cols(currentData,labels,'spos','cell num',cellNum,'session',sessions{1},'exptHalf',expHalves{1});

    subplot(1,2,1)
    plot(posAve(:,2),posAve(:,3),'-','color',[0.5 0.5 0.5])
    hold on
    plot(spos(:,2),spos(:,3),'k.','markersize',15)
    set(gca,'ydir','reverse')
    axis off
    axis equal
    
    if i == 3
        iSession = 4;
    else
        iSession = 2;
    end
    
    posAve = extract.cols(currentData,labels,'p','cell num',cellNum,'session',sessions{iSession},'exptHalf',expHalves{1});
    spikes = extract.cols(currentData,labels,'s','cell num',cellNum,'session',sessions{iSession},'exptHalf',expHalves{1});
    spos = extract.cols(currentData,labels,'spos','cell num',cellNum,'session',sessions{iSession},'exptHalf',expHalves{1});
    
    subplot(1,2,2)
    plot(posAve(:,2),posAve(:,3),'-','color',[0.5 0.5 0.5])
    hold on
    plot(spos(:,2),spos(:,3),'k.','markersize',15)
    set(gca,'ydir','reverse')
    axis off
    axis equal
    
    figure;
    mapBL1 = extract.cols(currentData,labels,'rate map','cell num',cellNum,'exptHalf',expHalves{1},'session',sessions{1});
    peakBL1 = extract.cols(currentData,labels,'peak rate','cell num',cellNum,'exptHalf',expHalves{1},'session',sessions{1});
    subplot(1,2,1)
    colorMapBRK(mapBL1);
    title(sprintf('Cell %s',cellNum))
    text(0,35,sprintf('%d Hz',round(peakBL1)))
    
    mapCNO1 = extract.cols(currentData,labels,'rate map','cell num',cellNum,'exptHalf',expHalves{1},'session',sessions{iSession});
    peakCNO1 = extract.cols(currentData,labels,'peak rate','cell num',cellNum,'exptHalf',expHalves{1},'session',sessions{iSession});
    subplot(1,2,2)
    colorMapBRK(mapCNO1);
    text(0,35,sprintf('%d Hz',round(peakCNO1)))
    
    sc(i,1) = analyses.spatialCrossCorrelation(mapBL1,mapCNO1);
end

%%
startup
CYL
hippoGlobe.arena
load masterMatAll_v4_2026

currentData = dataOutput;
sessions = unique(extract.cols(currentData,labels,'session'),'stable');
clusterNums = unique(extract.cols(currentData,labels,'cell num'),'stable');

%%
cellNum = 67;
cellNum = num2str(cellNum);

mapBL1 = extract.cols(currentData,labels,'rate map','cell num',cellNum,'session','A1');
peakBL1 = extract.cols(currentData,labels,'peak rate','cell num',cellNum,'session','A1');

figure;
subplot(1,2,1)
colorMapBRK(mapBL1);
title(sprintf('Cell %s',cellNum))
text(0,35,sprintf('%d Hz',round(peakBL1)))

mapCNO1 = extract.cols(currentData,labels,'rate map','cell num',cellNum,'session','B1');
peakCNO1 = extract.cols(currentData,labels,'peak rate','cell num',cellNum,'session','B1');
subplot(1,2,2)
colorMapBRK(mapCNO1);
text(0,35,sprintf('%d Hz',round(peakCNO1)))

scNov = analyses.spatialCrossCorrelation(mapBL1,mapCNO1);

posAve = extract.cols(currentData,labels,'p','cell num',cellNum,'session',sessions{1});
spikes = extract.cols(currentData,labels,'s','cell num',cellNum,'session',sessions{1});
spos = extract.cols(currentData,labels,'spos','cell num',cellNum,'session',sessions{1});

figure;
subplot(1,2,1)
plot(posAve(:,2),posAve(:,3),'-','color',[0.5 0.5 0.5])
hold on
plot(spos(:,2),spos(:,3),'k.','markersize',15)
set(gca,'ydir','reverse')
axis off
axis equal

posAve = extract.cols(currentData,labels,'p','cell num',cellNum,'session',sessions{2});
spikes = extract.cols(currentData,labels,'s','cell num',cellNum,'session',sessions{2});
spos = extract.cols(currentData,labels,'spos','cell num',cellNum,'session',sessions{2});

subplot(1,2,2)
plot(posAve(:,2),posAve(:,3),'-','color',[0.5 0.5 0.5])
hold on
plot(spos(:,2),spos(:,3),'k.','markersize',15)
set(gca,'ydir','reverse')
axis off
axis equal
