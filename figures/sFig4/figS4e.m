load figS4e

% BL1xBL2
figure;
[f,x,flo,fup] = ecdf(conAll(:,2));
[h1,h2,h3] = shadedplot(x',flo',fup',[252/255 186/255 131/255]);
set(h1,'FaceAlpha',0.3,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(conAll(:,2));
set(C,'color',[252/255 186/255 131/255],'linewidth',3)

hold on
[f,x,flo,fup] = ecdf(exp1x(:,2));
[h1,h2,h3] = shadedplot(x',flo',fup',[166/255 189/255 220/255]);
set(h1,'FaceAlpha',0.3,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(exp1x(:,2));
set(C,'color',[166/255 189/255 220/255],'linewidth',3)

xlabel('Spatial correlation')
ylabel('Cumulative proprotion')
title ''
grid off
set(gca,'ytick',0:0.25:1)

[~,p,ks] = kstest2(exp1x(:,2),conAll(:,2))
sum(~isnan(exp1x(:,2)))
sum(~isnan(conAll(:,2)))

nanmedian(exp1x(:,2))
nanmedian(conAll(:,2))