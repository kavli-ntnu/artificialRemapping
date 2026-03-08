load fig1b

xVals = [1 1.5 2.25 2.75 3.5 4];
figure;
plotSpread({conDiff(:,1) expDiff(:,1) conDiff(:,2) expDiff(:,2) conDiff(:,3) expDiff(:,3)},'xvalues',xVals);
o = flipud(findobj(gca,'type','line'));
set(o([1 3 5]),'markersize',8,'color',[0.8 0.8 0.8])
set(o([2 4 6]),'markersize',8,'color',[0.5 0.5 0.5])

set(gca,'xtickmode','manual','xtick',[1.25 2.5 3.75],'xticklabels',{'BL-BL','BL-CNO 1','BL-CNO 2'},'box','off','fontsize',14,'fontweight','bold')
ylabel('Subfield rate change','fontsize',14,'fontweight','bold')

nanmedian(conDiff)
nanmedian(expDiff)

sum(~isnan(conDiff))
sum(~isnan(expDiff))

[p,~,s] = ranksum(conDiff(:,1),expDiff(:,1))
[p,~,s] = ranksum(conDiff(:,2),expDiff(:,2),'tail','left')
[p,~,s] = ranksum(conDiff(:,3),expDiff(:,3),'tail','left')