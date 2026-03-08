load figS4b_4d

figure;
set(gcf,'color','w')
hold on
xVals = [0.5 1.25 2.5 3.25];

h = bar(xVals(1:2:end),nanmedian(conAllAbs(:,2:3)));
set(h,'barwidth',0.25,'facecolor',[1 1 1],'edgecolor',[.8 .8 .8],'LineWidth',2.5)

h = bar(xVals(2:2:end),nanmedian(exp1xAbs(:,2:3)));
set(h,'barwidth',0.25,'facecolor',[1 1 1],'edgecolor',[.5 .5 .5],'LineWidth',2.5)

plotSpread({conAllAbs(:,2) exp1xAbs(:,2) conAllAbs(:,3) exp1xAbs(:,3)},'xvalues',xVals);
o = flipud(findobj(gca,'type','line'));
set(o([1 3]),'markersize',8,'color',[0.8 0.8 0.8])
set(o([2 4]),'markersize',8,'color',[0.5 0.5 0.5])

set(gca,'xtickmode','manual','xtick',[0.87 2.87],'xticklabels',{'BL-CNO 1','BL-CNO 2'},'box','off','fontsize',14,'fontweight','bold')
ylabel('Size change','fontsize',14,'fontweight','bold')
xlim([-0.25 4])

[p,~,stats] = ranksum(conAllAbs(:,2),exp1xAbs(:,2),'tail','left')
sum(~isnan(exp1xAbs(:,2)))
sum(~isnan(conAllAbs(:,2)))

[p,~,stats] = ranksum(conAllAbs(:,3),exp1xAbs(:,3),'tail','left')
sum(~isnan(exp1xAbs(:,3)))
sum(~isnan(conAllAbs(:,3)))

%% scatter plot of difference scores 
xVals = exp1x(:,2);
yVals = exp1x(:,3);
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

figure
set(gcf,'color','w')
plot(xVals,yVals,'k.','markersize',15)
hold on
axis([-1 1 -1 1])

[r,p] = corr(xVals,yVals)
r2 = calc.benFit(xVals,yVals,'degree',1,'showLine',1)

set(gca,'box','off','fontsize',14,'fontweight','bold')
xlabel('Size change BL-CNO 1','fontsize',14,'fontweight','bold')
ylabel('Size change BL-CNO 2','fontsize',14,'fontweight','bold')
text(-0.95,0.98,sprintf('r = %.2f',r),'fontsize',14,'fontweight','bold')
