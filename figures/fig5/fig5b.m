load fig5b_5c

%%
cnt = 1;
figure;
selectedCells = [28 108 123];
for i = 1:size(selectedCells,2)
    cellNum = selectedCells(i);
    
    mapBL = pBL(:,:,cellNum);
    subplot(3,3,cnt)
    colorMapBRK(general.smooth(mapBL,6.6667));
    title(sprintf('%d',cellNum))
    cnt = cnt + 1;
    
    mapCNO1 = pCNO1(:,:,cellNum);
    subplot(3,3,cnt)
    colorMapBRK(general.smooth(mapCNO1,6.6667));
    cnt = cnt + 1;
    
    mapCNO2 = pCNO2(:,:,cellNum);
    subplot(3,3,cnt)
    colorMapBRK(general.smooth(mapCNO2,6.6667));
    cnt = cnt + 1;
    
    sc(i,1) = analyses.spatialCrossCorrelation(general.smooth(mapCNO1,6.6667),general.smooth(mapCNO2,6.6667));
end

