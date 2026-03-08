load fig2a

figure('position',[614,-3,560,866]);
cnt = 1;
for i = 1:size(mapStore,1)
    for j = 1:size(mapStore,2)
        subplot(6,size(mapStore,2),cnt)
        
        % scale to peaks when there are large rate changes, otherwise not
        if i == 2 && j == 1
            colorMapBRK(mapStore{i,j},'cutoffs',[0 peakRates(i,2)]);
        elseif i == 2 && j == 3
            colorMapBRK(mapStore{i,j},'cutoffs',[0 peakRates(i,4)]);
        elseif i == 3 && j == 2
            colorMapBRK(mapStore{i,j},'cutoffs',[0 peakRates(i,1)]);
        elseif i == 3 && j == 4
            colorMapBRK(mapStore{i,j},'cutoffs',[0 peakRates(i,3)]);
        else
            colorMapBRK(mapStore{i,j});
        end
        
        text(0,size(mapStore{i,j},1)+3,sprintf('%d Hz',round(peakRates(i,j))))
        cnt = cnt + 1;
    end
    
    sc(i,1) = analyses.spatialCrossCorrelation(mapStore{i,2},mapStore{i,4});
end

