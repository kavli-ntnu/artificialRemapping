load figS6d

figure;
hold on

[f,x,flo,fup] = ecdf(sc_realign);
shadedplot(x',flo',fup',[0.5 0.5 0.5],'k');
C = cdfplot(sc_realign);
set(C,'color','k','linewidth',0.5)

[f,x,flo,fup] = ecdf(sc_subfieldRate);
shadedplot(x',flo',fup','w','k');
C = cdfplot(sc_subfieldRate);
set(C,'color','k','linewidth',0.5)

title ''
grid off
box off
ylabel 'Proportion'
xlabel 'Spatial correlation'
set(gca,'fontsize',14,'fontw','bold')
xlim([-0.6 1])

sum(~isnan(sc_realign))
sum(~isnan(sc_subfieldRate))

nanmedian(sc_realign)
nanmedian(sc_subfieldRate)

[~,p,k] = kstest2(sc_realign,sc_subfieldRate)

