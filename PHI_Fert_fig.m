% Create figure illustrating different fertilization curves, varying the
% Phi (male importance) parameter

m = 0:0.01:1; % sex ratio
a = 1 ; % Beta parameter

figure
hold on

pf = betacdf(m,a,2);
ph1 = plot(m,pf);

pf = betacdf(m,a,6);
ph2 = plot(m,pf);

pf = betacdf(m,a,12);
ph3 = plot(m,pf);

pf = betacdf(m,a,20);
ph4 = plot(m,pf);

set(ph4, 'color','k','linewidth',1.5)
set(ph3,'color',[0.8 0.2 0.2],'linewidth',1.5)
set(ph2,'color',[0.7 0.2 0.5],'linewidth',1.5)
set(ph1,'color',[0.2 0.2 0.8],'linewidth',1.5)

set(gca,'tickdir','out','ticklength',[0.02 0.02]) 
set(gca,'linewidth',1,'fontsize',14)
set(gca,'xtick',0:0.1:1)
set(gca,'ytick',0:0.1:1)
set(gca,'ylim',[0 1.0])
set(gca,'xlim',[0 1.0]) 

xlabel(gca,'Sex ratio (proportional male biomass)','fontsize',18)
ylabel(gca,'Fertilization rate','fontsize',18)
%title('Network-Persistence')
legend([ph1,ph2,ph3,ph4],'Phi = 2','Phi = 6','Phi = 12','Phi = 20','Location','northeast') ;

