%%
load fig4c

figure;
binRange = 0:3.5:100;
centers = binRange(1:end-1) + diff(binRange)/2;
counts = histcounts(errorVals,binRange) / sum(~isnan(errorVals));
bar(centers,counts,'barwidth',1)
xticks(20:20:100)
% set(gca,'xtickmode','manual','xtick',0:10:100,'box','off','fontsize',14,'fontweight','bold')
xlabel('Prediction offset (cm)')
ylabel('Proportion of neurons')
ylim([0 0.15])

sum(~isnan(errorVals))
nanmedian(errorVals)

%%
load fig4d

figure;
binRange = 0:3.5:100;
centers = binRange(1:end-1) + diff(binRange)/2;
counts = histcounts(shuffleVals,binRange) / sum(~isnan(shuffleVals));
bar(centers,counts,'barwidth',1)
xticks(20:20:100)
xlabel('Prediction offset (cm)')
ylabel('Proportion of neurons')
ylim([0 0.15])

sum(~isnan(shuffleVals))
nanmedian(shuffleVals)

[p,~,stats] = ranksum(errorVals,shuffleVals,'tail','left')

%%
[f,xi] = ksdensity(errorVals,0:1:100);
figure;
plot(xi,f,'color','b')
hold on

[f,xi] = ksdensity(shuffleVals,0:1:100);
plot(xi,f,'color','g')

xlabel('Prediction offset (cm)')
ylabel('Probability density')
box off
ylim([0 0.0235])
set(gca,'xtick',0:20:100)
