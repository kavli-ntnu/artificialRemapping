%%
load fig3b_3c_3d

% hM3
binRange = 0:2:60;
counts = histcounts(err,binRange) / sum(~isnan(err));
figure;
bar(1:2:59,counts,'barwidth',1,'edgecolor','k','facecolor',[0.5 0.5 0.5])
set(gca,'xtickmode','manual','xtick',0:10:60,'box','off','fontsize',14,'fontweight','bold')
xlabel('Prediction offset (cm)')
ylabel('Proportion of neurons')
ylim([0 0.16])

nanmedian(err)
sum(~isnan(err))

%% AxB
binRange = 0:2:60;
counts = histcounts(nov,binRange) / sum(~isnan(nov));
figure;
bar(1:2:59,counts,'barwidth',1,'edgecolor','k','facecolor',[0.5 0.5 0.5])
set(gca,'xtickmode','manual','xtick',0:10:60,'box','off','fontsize',14,'fontweight','bold')
xlabel('Prediction offset (cm)')
ylabel('Proportion of neurons')
ylim([0 0.16])

nanmedian(nov)
sum(~isnan(nov))

%%
[p,~,stats] = ranksum(err,nov,'tail','left')

%%
nBins = 50;
figure;
h1 = histfit(err,nBins,'kernel');
x1 = h1(2).XData;
y1 = h1(2).YData;
counts1 = histcounts(err,nBins);
close(gcf)

figure;
h2 = histfit(nov,nBins,'kernel');
x2 = h2(2).XData;
y2 = h2(2).YData;
counts2 = histcounts(nov,nBins);
close(gcf)

figure;
plot(x1,y1/sum(counts1),'b-')
hold on
plot(x2,y2/sum(counts2),'g-')
xlabel('Prediction offset (cm)')
ylabel('Probability density')
ylim([0 0.042])
xlim([-0.02 60])
set(gca,'xtick',0:10:60,'ytick',0:0.01:0.05)




