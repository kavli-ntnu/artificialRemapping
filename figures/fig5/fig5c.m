load fig5b_5c

%% unsmoothed
spatialCorrAll = nan(200*numSims,3);
for i = 1:200*numSims
    if iBL(i) || iCNO1(i)
        mapBL = pBL(:,:,i);
        mapCNO = pCNO1(:,:,i);
        spatialCorrAll(i,1) = analyses.spatialCrossCorrelation(mapBL,mapCNO);
    end
    
    if iBL(i) || iCNO2(i)
        mapBL = pBL(:,:,i);
        mapCNO = pCNO2(:,:,i);
        spatialCorrAll(i,2) = analyses.spatialCrossCorrelation(mapBL,mapCNO);
    end
    
    if iCNO1(i) || iCNO2(i)
        mapCNO1 = pCNO1(:,:,i);
        mapCNO2 = pCNO2(:,:,i);
        spatialCorrAll(i,3) = analyses.spatialCrossCorrelation(mapCNO1,mapCNO2);
    end 
end

%%
figure;
[f,x,flo,fup] = ecdf(spatialCorrAll(:,1));
[h1,h2,h3] = shadedplot(x',flo',fup',[10/255 113/255 178/255]);
set(h1,'FaceAlpha',0.2,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(spatialCorrAll(:,1));
set(C,'color',[10/255 113/255 178/255 0.4],'linewidth',3)

hold on
[f,x,flo,fup] = ecdf(spatialCorrAll(:,2));
[h1,h2,h3] = shadedplot(x',flo',fup',[10/255 113/255 178/255]);
set(h1,'FaceAlpha',0.3,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(spatialCorrAll(:,2));
set(C,'color',[10/255 113/255 178/255],'linewidth',3)

hold on
[f,x,flo,fup] = ecdf(spatialCorrAll(:,3));
[h1,h2,h3] = shadedplot(x',flo',fup',[0.7 0.7 0.7]);
set(h1,'FaceAlpha',0.3,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(spatialCorrAll(:,3));
set(C,'color','k','linewidth',3)

xlabel('Spatial correlation')
ylabel('Proprotion')
title ''
grid off
set(gca,'ytick',0:0.25:1)
xlim([-0.2 1])

%% smoothed
spatialCorrAll = nan(200*numSims,3);
for i = 1:200*numSims
    if iBL(i) || iCNO1(i)
        mapBL = general.smooth(pBL(:,:,i),6.67);
        mapCNO = general.smooth(pCNO1(:,:,i),6.67);
        spatialCorrAll(i,1) = analyses.spatialCrossCorrelation(mapBL,mapCNO);
    end
    
    if iBL(i) || iCNO2(i)
        mapBL = general.smooth(pBL(:,:,i),6.67);
        mapCNO = general.smooth(pCNO2(:,:,i),6.67);
        spatialCorrAll(i,2) = analyses.spatialCrossCorrelation(mapBL,mapCNO);
    end
    
    if iCNO1(i) || iCNO2(i)
        mapCNO1 = general.smooth(pCNO1(:,:,i),6.67);
        mapCNO2 = general.smooth(pCNO2(:,:,i),6.67);
        spatialCorrAll(i,3) = analyses.spatialCrossCorrelation(mapCNO1,mapCNO2);
    end 
end

nanmedian(spatialCorrAll)
sum(~isnan(spatialCorrAll))

%%
figure;
[f,x,flo,fup] = ecdf(spatialCorrAll(:,1));
[h1,h2,h3] = shadedplot(x',flo',fup',[10/255 113/255 178/255]);
set(h1,'FaceAlpha',0.2,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(spatialCorrAll(:,1));
set(C,'color',[10/255 113/255 178/255 0.4],'linewidth',3)

hold on
[f,x,flo,fup] = ecdf(spatialCorrAll(:,2));
[h1,h2,h3] = shadedplot(x',flo',fup',[10/255 113/255 178/255]);
set(h1,'FaceAlpha',0.3,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(spatialCorrAll(:,2));
set(C,'color',[10/255 113/255 178/255],'linewidth',3)

hold on
[f,x,flo,fup] = ecdf(spatialCorrAll(:,3));
[h1,h2,h3] = shadedplot(x',flo',fup',[0.7 0.7 0.7]);
set(h1,'FaceAlpha',0.3,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(spatialCorrAll(:,3));
set(C,'color','k','linewidth',3)

xlabel('Spatial correlation')
ylabel('Proprotion')
title ''
grid off
set(gca,'ytick',0:0.25:1)
xlim([-0.2 1])

%%
sum(~isnan(spatialCorrAll))

[~,p,ks] = kstest2(spatialCorrAll(:,1),spatialCorrAll(:,2))
[~,p,ks] = kstest2(spatialCorrAll(:,1),spatialCorrAll(:,3))
[~,p,ks] = kstest2(spatialCorrAll(:,2),spatialCorrAll(:,3))
