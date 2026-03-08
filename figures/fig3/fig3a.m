load fig3a

figure('position',[614,-3,560,866]);
cnt = 1;
for i = 1:size(mapStore,1)
    for j = 1:size(mapStore,2)
        subplot(size(mapStore,1),size(mapStore,2),cnt)
        colorMapBRK(mapStore{i,j});
        text(0,size(mapStore{i,j},1)+3,sprintf('%d Hz',round(peakRates(i,j))))
        cnt = cnt + 1;
    end
    
    sc(i,1) = analyses.spatialCrossCorrelation(mapStore{i,1},mapStore{i,2});
end