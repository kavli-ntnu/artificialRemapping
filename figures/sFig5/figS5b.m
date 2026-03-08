load figS5b

binRange = 0:2:60;
counts = histcounts(shuffleVals,binRange) / sum(~isnan(shuffleVals));
figure;
bar(1:2:59,counts,'barwidth',1,'edgecolor','k','facecolor',[0.5 0.5 0.5])
set(gca,'xtickmode','manual','xtick',0:10:60,'box','off','fontsize',14,'fontweight','bold')
xlabel('Prediction offset (cm)')
ylabel('Proportion of neurons')
ylim([0 0.16])

sum(~isnan(shuffleVals))
nanmedian(shuffleVals)

load fig3b_3c_3d
[p,~,stats] = ranksum(err,shuffleVals,'tail','left')
