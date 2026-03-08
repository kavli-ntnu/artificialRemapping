load figS5c

figure;
plotSpread({rdCon rdExp})
set(gca,'xtick',1:2,'xticklabels',{'Con','hM3'})
set(gca,'ytick',-1:0.5:1)

title('Rate change within BL primary field')

nanmedian(rdCon)
nanmedian(rdExp)

sum(~isnan(rdCon))
sum(~isnan(rdExp))

[p,~,s] = ranksum(rdCon,rdExp)


%%
load figS5d

figure;
plotSpread({rdCon rdExp})
set(gca,'xtick',1:2,'xticklabels',{'Con','hM3'})
set(gca,'ytick',-1:0.5:1)

title('Rate change within CNO primary field')

nanmedian(rdCon)
nanmedian(rdExp)

sum(~isnan(rdCon))
sum(~isnan(rdExp))

[p,~,s] = ranksum(rdCon,rdExp)