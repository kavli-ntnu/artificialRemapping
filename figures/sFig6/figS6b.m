load figS6b

figure;
hold on

[f,x,flo,fup] = ecdf(hM3);
shadedplot(x',flo',fup',[102/255 178/255 1],[0 0 1]);
C = cdfplot(hM3);
set(C,'color','b','linewidth',2)

[f,x,flo,fup] = ecdf(sc6);
shadedplot(x',flo',fup',[1 .5 .5],[1 0 0]);
C = cdfplot(sc6);
set(C,'color','r','linewidth',2)

[f,x,flo,fup] = ecdf(sc3);
shadedplot(x',flo',fup',[204/255 1 229/255],[0 153/255 76/255]);
C = cdfplot(sc3);
set(C,'color',[0 153/255 76/255],'linewidth',2)

[f,x,flo,fup] = ecdf(sc9);
shadedplot(x',flo',fup',[1 229/255 204/255],[1 128/255 0]);
C = cdfplot(sc9);
set(C,'color',[1 128/255 0],'linewidth',2)

title ''
grid off
box off
ylabel 'Cumulative proportion'
xlabel 'Spatial correlation'
set(gca,'fontsize',14,'fontw','bold','ytick',[0 .25 .5 .75 1])
xlim([-0.5 1])

sum(~isnan(hM3))
sum(~isnan(sc3))
sum(~isnan(sc6))
sum(~isnan(sc9))

[p,~,s] = ranksum(hM3,sc6)

nanmedian(sc3)
nanmedian(sc6)
nanmedian(sc9)

