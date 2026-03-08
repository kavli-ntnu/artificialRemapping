load figS3b_3d

% bar graph + plot spread
figure;
set(gcf,'color','w')
hold on
xVals = [0.5 1.25 2.5 3.25];

h = bar(xVals(1:2:end),nanmedian(conAllAbs(:,1:2)));
% set(h,'barwidth',0.3,'facecolor',[.8 .8 .8])
set(h,'barwidth',0.25,'facecolor',[1 1 1],'edgecolor',[.8 .8 .8],'LineWidth',2.5)

h = bar(xVals(2:2:end),nanmedian(expAllAbs(:,1:2)));
% set(h,'barwidth',0.3,'facecolor',[.5 .5 .5])
set(h,'barwidth',0.25,'facecolor',[1 1 1],'edgecolor',[.5 .5 .5],'LineWidth',2.5)

plotSpread({conAllAbs(:,1) expAllAbs(:,1) conAllAbs(:,2) expAllAbs(:,2)},'xvalues',xVals);
o = flipud(findobj(gca,'type','line'));
set(o([1 3]),'markersize',8,'color',[0.8 0.8 0.8])
set(o([2 4]),'markersize',8,'color',[0.5 0.5 0.5])

set(gca,'xtickmode','manual','xtick',[0.87 2.87],'xticklabels',{'BL-CNO 1','BL-CNO 2'},'box','off','fontsize',14,'fontweight','bold')
ylabel('Size change','fontsize',14,'fontweight','bold')
xlim([-0.25 4])

[p,~,stats] = ranksum(expAllAbs(:,1),conAllAbs(:,1),'tail','right')
sum(~isnan(conAllAbs(:,1)))
sum(~isnan(expAllAbs(:,1)))

[p,~,stats] = ranksum(expAllAbs(:,2),conAllAbs(:,2),'tail','right')
sum(~isnan(conAllAbs(:,2)))
sum(~isnan(expAllAbs(:,2)))

%% scatter plot of rate difference scores
xVals = expAll(:,1);
yVals = expAll(:,2);
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

figure
set(gcf,'color','w')
plot(xVals,yVals,'k.','markersize',15)
hold on

[r,p] = corr(xVals,yVals)
r2 = calc.benFit(xVals,yVals,'degree',1,'showLine',1)
sum(~isnan(xVals))

set(gca,'box','off','fontsize',14,'fontweight','bold')
xlabel('Size change Day 1','fontsize',14,'fontweight','bold')
ylabel('Size change Day 2','fontsize',14,'fontweight','bold')

text(-0.45,0.98,sprintf('r = %.2f',r),'fontsize',14,'fontweight','bold')
