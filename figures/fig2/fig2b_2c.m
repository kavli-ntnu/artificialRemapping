load fig2b_2c

% BLxCNO1 & BLxCNO2
figure;
[f,x,flo,fup] = ecdf(conAll(:,3));
[h1,h2,h3] = shadedplot(x',flo',fup',[219/255 61/255 49/255]);
set(h1,'FaceAlpha',0.2,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(conAll(:,3));
set(C,'color',[219/255 61/255 49/255 0.4],'linewidth',3)

hold on
[f,x,flo,fup] = ecdf(conAll(:,4));
[h1,h2,h3] = shadedplot(x',flo',fup',[219/255 61/255 49/255]);
set(h1,'FaceAlpha',0.3,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(conAll(:,4));
set(C,'color',[219/255 61/255 49/255],'linewidth',3)

[f,x,flo,fup] = ecdf(exp1x(:,3));
[h1,h2,h3] = shadedplot(x',flo',fup',[10/255 113/255 178/255]);
set(h1,'FaceAlpha',0.2,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(exp1x(:,3));
set(C,'color',[10/255 113/255 178/255 0.4],'linewidth',3)

hold on
[f,x,flo,fup] = ecdf(exp1x(:,4));
[h1,h2,h3] = shadedplot(x',flo',fup',[10/255 113/255 178/255]);
set(h1,'FaceAlpha',0.3,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(exp1x(:,4));
set(C,'color',[10/255 113/255 178/255],'linewidth',3)

xlabel('Spatial correlation')
ylabel('Cumulative proprotion')
title ''
grid off
set(gca,'ytick',0:0.25:1)

%% CNO1xCNO2
figure;
[f,x,flo,fup] = ecdf(conAll(:,5));
[h1,h2,h3] = shadedplot(x',flo',fup',[132/255 32/255 34/255]);
set(h1,'FaceAlpha',0.3,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(conAll(:,5));
set(C,'color',[132/255 32/255 34/255],'linewidth',3)

hold on
[f,x,flo,fup] = ecdf(exp1x(:,5));
[h1,h2,h3] = shadedplot(x',flo',fup',[25/255 72/255 101/255]);
set(h1,'FaceAlpha',0.3,'EdgeColor','none')
set(h2,'color','none')
set(h3,'color','none')
C = cdfplot(exp1x(:,5));
set(C,'color',[25/255 72/255 101/255],'linewidth',3)

xlabel('Spatial correlation')
ylabel('Cumulative proprotion')
title ''
grid off
set(gca,'ytick',0:0.2:1)

%% BL1xBL2
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

%%
sum(~isnan(exp1x))
sum(~isnan(conAll))

[~,p,ks] = kstest2(exp1x(:,3),conAll(:,3))
[~,p,ks] = kstest2(exp1x(:,4),conAll(:,4))
[~,p,ks] = kstest2(exp1x(:,5),conAll(:,5))
