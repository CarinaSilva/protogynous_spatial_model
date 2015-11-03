function Network_MinCR_persist_fig(Lf,f)

% Plot persistence thresholds for network persistence scenario

% Inputs: Lf = average size at entry to fishery
%         f = 1 or 2, corresponding to FLEP = [0, 0.2] in simulations

load(strcat('spatialnetwork_minRS_Lf',num2str(Lf),'_persist.mat'),'PHI','MinRS')


figure
clf
hold on

ph1 = plot(PHI,squeeze(MinRS(1,:,f))) ;
ph2 = plot(PHI,squeeze(MinRS(2,:,f))) ;
ph3 = plot(PHI,squeeze(MinRS(3,:,f))) ;
ph4 = plot(PHI,squeeze(MinRS(4,:,f))) ;

set(ph1, 'color','k','linewidth',1.5)
set(ph2,'color',[0.8 0.2 0.2],'linewidth',1.5)
set(ph3,'color',[0.7 0.2 0.5],'linewidth',1.5)
set(ph4,'color',[0.2 0.2 0.8],'linewidth',1.5)


set(gca,'tickdir','out','ticklength',[0.02 0.02]) % gca = 'get current axis'
set(gca,'linewidth',1,'fontsize',14)
set(gca,'xtick',1:2:20)
set(gca,'ytick',0:0.1:1.0)
xlabel(gca,'High              Male importance (   )               Low','fontsize',16)
ylabel(gca,{'Minimum fraction of coastline'; 'required for persistence'},'fontsize',16)
%legend([ph1,ph2,ph3,ph4],'GON','r = 1','r = 10','r = 100','Location','northeast')
set(gca,'ylim',[0 0.7])
set(gca,'xlim',[1 20])
%title('Network-Persistence')
legend([ph1,ph2,ph3,ph4],'Gonochore','Fixed','Relative Size','Size Frequency','Location','northeast');