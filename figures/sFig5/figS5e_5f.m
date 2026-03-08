% to load source data:
load figS5e_5f

% to generate source data, use: fx_fig3b_3c_3d

%% scatter plot of prediction error vs firing rate in predicted subfield
xVals = err;
yVals = peakRate;
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

figure
set(gcf,'color','w')
scatter(xVals,yVals,'k','lineWidth',1.5)
hold on

[r_val p_val] = corr(xVals,yVals)

set(gca,'box','off','fontsize',14,'fontweight','bold')
xlabel('Prediction offset (cm)','fontsize',14,'fontweight','bold')
ylabel('Peak rate (Hz)','fontsize',14,'fontweight','bold')
text(40,13.5,sprintf('r = %.2f',r_val),'fontsize',14,'fontweight','bold')
text(40,12.7,sprintf('p = %.2f',p_val),'fontsize',14,'fontweight','bold')

%% scatter plot of prediction error vs spatial corr
xVals = err;
yVals = spatialCorr;
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

figure
set(gcf,'color','w')
scatter(xVals,yVals,'k','lineWidth',1.5)
hold on

[r_val p_val] = corr(xVals,yVals)

set(gca,'box','off','fontsize',14,'fontweight','bold')
xlabel('Prediction offset','fontsize',14,'fontweight','bold')
ylabel('Spatial correlation BL-CNO','fontsize',14,'fontweight','bold')
text(19,1,sprintf('r = %.2f',r_val),'fontsize',14,'fontweight','bold')

