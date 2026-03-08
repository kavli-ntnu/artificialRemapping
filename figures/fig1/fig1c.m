load fig1c

x = reshape(diffDay1,[],1);
y = reshape(diffDay2,[],1);
xnan = isnan(x);
ynan = isnan(y);
x = x(~xnan & ~ynan);
y = y(~xnan & ~ynan);

figure;
set(gcf,'color','w')
scatter(x,y,30,[0 0 0],'filled')
hold on

r2 = calc.benFit(x,y,'degree',1,'showLine',1)
[r,p] = corr(x,y)
sum(~isnan(x))

set(gca,'box','off','fontsize',14,'fontweight','bold')
xlabel('Subfield rate change BL-CNO 1','fontsize',14,'fontweight','bold')
ylabel('Subfield rate change BL-CNO2','fontsize',14,'fontweight','bold')
xlim([-1 1])
ylim([-1 1])
text(-0.95,0.95,sprintf('hM3: r = %.2f, p = %.1e',r, p),'color','k','fontweight','bold','fontsize',14) 

%%
figure;
expAll = randAll;

ratesBL1 = reshape(expAll(:,:,1),[],1);
ratesCNO1 = reshape(expAll(:,:,2),[],1);

diffDay1 = (ratesCNO1 - ratesBL1) ./ (ratesCNO1 + ratesBL1);

ratesBL2 = reshape(expAll(:,:,3),[],1);
ratesCNO2 = reshape(expAll(:,:,4),[],1);

diffDay2 = (ratesCNO2 - ratesBL2) ./ (ratesCNO2 + ratesBL2);
x = reshape(diffDay1,[],1);
y = reshape(diffDay2,[],1);
xnan = isnan(x);
ynan = isnan(y);
x = x(~xnan & ~ynan);
y = y(~xnan & ~ynan);

set(gcf,'color','w')
scatter(x,y,30,[189/255 189/255 189/255],'filled')
hold on

r2 = calc.benFit(x,y,'degree',1,'showLine',1) %Need to go in the function to change the line color
[r,p] = corr(x,y)
sum(~isnan(x))

set(gca,'box','off','fontsize',14,'fontweight','bold')
xlabel(' grid subfield rate BL-CNO 1','fontsize',14,'fontweight','bold')
ylabel(' grid subfield rate BL-CNO2','fontsize',14,'fontweight','bold')
% axis square
xlim([-1 1])
ylim([-1 1])

text(-0.95,0.85,sprintf('Shuffled: r = %.2f, p = %.1e',r, p),'color','#969696','fontweight','bold','fontsize',14)

