load figS3e

figure;
set(gcf,'color','w')
hold on
xVals = [0.5 1.25 2.5 3.25];

plotSpread({conAll(:,1) expAll(:,1) conAll(:,2) expAll(:,2)},'xvalues',xVals);
o = flipud(findobj(gca,'type','line'));
set(o([1 3]),'markersize',8,'color',[0.8 0.8 0.8])
set(o([2 4]),'markersize',8,'color',[0.5 0.5 0.5])

set(gca,'xtickmode','manual','xtick',[0.87 2.87],'xticklabels',{'BL-CNO 1','BL-CNO 2'},'box','off','fontsize',14,'fontweight','bold')
ylabel('Spatial correlation','fontsize',14,'fontweight','bold')
xlim([-0.25 4])

%%
[p,~,stats] = ranksum(conAll(:,1),expAll(:,1))
sum(~isnan(conAll(:,1)))
sum(~isnan(expAll(:,1)))

%%
[p,~,stats] = ranksum(conAll(:,2),expAll(:,2))
sum(~isnan(conAll(:,2)))
sum(~isnan(expAll(:,2)))