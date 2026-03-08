load figS6a

figure;
plotSpread({conBL_CNO expBL_CNO})
hold on
plot(1:2,[nanmedian(conBL_CNO) nanmedian(expBL_CNO)],'r+')
ylabel('Grid subfield rate changes')
set(gca,'ytick',0:0.2:1)
ylim([-0.02 1.02])

nanmedian(expBL_CNO)

sum(~isnan(conBL_CNO))
sum(~isnan(expBL_CNO))

[p,~,s] = ranksum(conBL_CNO,expBL_CNO)
