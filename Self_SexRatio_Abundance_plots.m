function Self_SexRatio_Abundance_plots

% Plot both Sex Ratio & Abundance for example scenarios

% Run the scenarios:
doRuns = false;
if doRuns    
    FLEPs = 0.3:-0.025:0;
    PHIs = [2 6 20];
    for f = 1:length(FLEPs)
        for p = 1:length(PHIs)
    Spatial_Struct(PHIs(p),FLEPs(f),'self');
        end
    end
end % end doRuns

Mid = 45; % where to split the coastline to plot the MPA in the middle

Col = [0 0 0; 0.8 0.2 0.2; 0.7 0.2 0.5; 0.2 0.2 0.8];
LW = 1.5;
Style = {'-','--','-.',':'};

%%%%%%%%%
% Figure 1: baseline case with Phi = 6
figure(1)
set(gcf,'units','cent','position',[10 10 16.9 15])
clf

% Panel B: sustainable fishing, sex ratio
load('spatialself_Oct2015_PHI6_FLEP1_baseline.mat','S','Spatial') ;
Spatial_U = Spatial;

load('spatialself_Oct2015_PHI6_FLEP0.2_baseline.mat','S','Spatial') ;
%------------------------------------------------------------------------
subplot(2,2,2);
hold on
%plot_patches

for s = 1:length(S)
Str = Spatial.(S{s}).F(1).PHI(1).RW(1).SexratioNum;
X(:,s) = [Str((Mid+1):end), Str(1:Mid)];
ph(s) = plot(X(:,s));
set(ph(s),'color',Col(s,:),'linewidth',LW,'linestyle',Style{s})
end

Ax = gca;
setup_axes(Ax)
ylabel(Ax,'Sex ratio (prop. male)','fontsize',10)
text(1,0.95,'B','fontsize',12)
%------------------------------------------------------------------------

%------------------------------------------------------------------------
subplot(2,2,1);
hold on
%plot_patches

for s = 1:length(S)
Str = Spatial.(S{s}).F(1).PHI(1).RW(1).BiomPatch;
Str = Str./Spatial_U.(S{s}).F(1).PHI(1).RW(1).BiomPatch; % scale by unfished
Str = Str(:)';
X(:,s) = [Str((Mid+1):end), Str(1:Mid)];
ph(s) = plot(X(:,s));
set(ph(s),'color',Col(s,:),'linewidth',LW,'linestyle',Style{s})
end

Ax = gca;
setup_axes(Ax)
ylabel(Ax,{'Male biomass'; '(relative to unfished)'},'fontsize',10)
text(1,0.95,'A','fontsize',12)
setup_legend(ph)
%------------------------------------------------------------------------

load('spatialself_march2015_PHI6_FLEP0.2_baseline.mat','S','Spatial') ;
%------------------------------------------------------------------------
subplot(2,2,4);
hold on
%plot_patches

for s = 1:length(S)
Str = Spatial.(S{s}).F(1).PHI(1).RW(1).SexratioNum;
X(:,s) = [Str((Mid+1):end), Str(1:Mid)];
ph(s) = plot(X(:,s));
set(ph(s),'color',Col(s,:),'linewidth',LW,'linestyle',Style{s})
end

Ax = gca;
setup_axes(Ax)
ylabel(Ax,'Sex ratio (prop. male)','fontsize',10)
text(1,0.95,'D','fontsize',12)
%------------------------------------------------------------------------

% Panel C: unsustainable fishing, Biomass
%------------------------------------------------------------------------
subplot(2,2,3);
hold on
%plot_patches

for s = 1:length(S)
Str = Spatial.(S{s}).F(1).PHI(1).RW(1).BiomPatch;
Str = Str./Spatial_U.(S{s}).F(1).PHI(1).RW(1).BiomPatch; % scale by unfished
Str = Str(:)';
X(:,s) = [Str((Mid+1):end), Str(1:Mid)];
ph(s) = plot(X(:,s));
set(ph(s),'color',Col(s,:),'linewidth',LW,'linestyle',Style{s})
end

Ax = gca;
setup_axes(Ax)
ylabel(Ax,{'Male biomass'; '(relative to unfished)'},'fontsize',10)
text(1,0.95,'C','fontsize',12)
%------------------------------------------------------------------------

%%%%%%%%%
% Figure 2: cases with Phi = 2, Phi = 20
figure(2)
set(gcf,'units','cent','position',[10 10 16.9 15])
clf

% Panel B: Phi = 2, unsustainable fishing, sex ratio
load('spatialself_march2015_PHI2_FLEP1_baseline.mat','S','Spatial') ;
Spatial_U = Spatial;

load('spatialself_march2015_PHI2_FLEP0.2_baseline.mat','S','Spatial') ;
%------------------------------------------------------------------------
subplot(2,2,2);
hold on
%plot_patches

for s = 1:length(S)
Str = Spatial.(S{s}).F(1).PHI(1).RW(1).SexratioNum;
X(:,s) = [Str((Mid+1):end), Str(1:Mid)];
ph(s) = plot(X(:,s));
set(ph(s),'color',Col(s,:),'linewidth',LW,'linestyle',Style{s})
end

Ax = gca;
setup_axes(Ax)
ylabel(Ax,'Sex ratio (prop. male)','fontsize',10)
text(1,0.95,'B','fontsize',12)

%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Panel A: Phi = 2, unsustainable fishing, population density
subplot(2,2,1);
hold on
%plot_patches

for s = 1:length(S)
Str = Spatial.(S{s}).F(1).PHI(1).RW(1).BiomPatch;
Str = Str./Spatial_U.(S{s}).F(1).PHI(1).RW(1).BiomPatch; % scale by unfished
Str = Str(:)';
X(:,s) = [Str((Mid+1):end), Str(1:Mid)];
ph(s) = plot(X(:,s));
set(ph(s),'color',Col(s,:),'linewidth',LW,'linestyle',Style{s})
end

Ax = gca;
setup_axes(Ax)
ylabel(Ax,{'Biomass'; '(relative to unfished)'},'fontsize',10)
text(1,0.95,'A','fontsize',12)
setup_legend(ph)
%------------------------------------------------------------------------

% Panel D: Phi = 20, unsustainable fishing, sex ratio
load('spatialself_march2015_PHI20_FLEP1_baseline.mat','S','Spatial') ;
Spatial_U = Spatial;
load('spatialself_march2015_PHI20_FLEP0.2_baseline.mat','S','Spatial') ;
%------------------------------------------------------------------------
subplot(2,2,4);
hold on
%plot_patches

for s = 1:length(S)
Str = Spatial.(S{s}).F(1).PHI(1).RW(1).SexratioNum;
X(:,s) = [Str((Mid+1):end), Str(1:Mid)];
ph(s) = plot(X(:,s));
set(ph(s),'color',Col(s,:),'linewidth',LW,'linestyle',Style{s})
end

Ax = gca;
setup_axes(Ax)
ylabel(Ax,'Sex ratio (prop. male)','fontsize',10)
text(1,0.95,'D','fontsize',12)

%------------------------------------------------------------------------

% Panel C: unsustainable fishing, Biomass
%------------------------------------------------------------------------
subplot(2,2,3);
hold on
%plot_patches

for s = 1:length(S)
Str = Spatial.(S{s}).F(1).PHI(1).RW(1).BiomPatch;
Str = Str./Spatial_U.(S{s}).F(1).PHI(1).RW(1).BiomPatch; % scale by unfished
Str = Str(:)';
X(:,s) = [Str((Mid+1):end), Str(1:Mid)];
ph(s) = plot(X(:,s));
set(ph(s),'color',Col(s,:),'linewidth',LW,'linestyle',Style{s})
end

Ax = gca;
setup_axes(Ax)
ylabel(Ax,{'Biomass'; '(relative to unfished)'},'fontsize',10)
text(1,0.95,'C','fontsize',12)
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Figure 3: Mean biomass & sex ratio vs. FLEP, for a range of values
figure(3)
set(gcf,'units','cent','position',[10 10 16.9 20])
clf
% 

PHIs = [20, 6, 2];

sPlots = [1, 2; 3, 4; 5, 6];

for p = 1:3
FLEPs = 0.3:-0.025:0;
load(strcat('spatialself_Oct2015_PHI',num2str(PHIs(p)),'_FLEP1_baseline.mat'),'S','Spatial') ;
Spatial_U = Spatial;
for f = 1:length(FLEPs)
load(strcat('spatialself_Oct2015_PHI',num2str(PHIs(p)),'_FLEP',num2str(FLEPs(f)),'_baseline.mat'),'S','Spatial') ;
Spatial_tmp(f) = Spatial;
end

% Plot Biomass
subplot(3,2,sPlots(p,1))
hold on
Str = [];
for s = 1:length(S)
    StrU = mean(Spatial_U.(S{s}).F(1).PHI(1).RS(1).BiomPatch);
    for f = 1:length(FLEPs)
    Str(f) = mean(Spatial_tmp(f).(S{s}).F(1).PHI(1).RS(1).BiomPatch);
    end 
Str_p = Str./StrU;
ph(s) = plot(FLEPs,Str_p);
set(ph(s),'color',Col(s,:),'linewidth',LW,'linestyle',Style{s})
end
set(gca,'tickdir','out','ticklength',[0.015 0.015])
ylim([0 1])

% Plot SexRatio
subplot(3,2,sPlots(p,2))
hold on
Str = [];
for s = 1:length(S)
    for f = 1:length(FLEPs)
    Str(f) = Spatial_tmp(f).(S{s}).F(1).PHI(1).RS(1).SexratioNum(end-5);
    end 
ph(s) = plot(FLEPs,Str);
set(ph(s),'color',Col(s,:),'linewidth',LW,'linestyle',Style{s})
end
set(gca,'tickdir','out','ticklength',[0.015 0.015])
ylim([0 1])
end % end loop over PHIs 
%------------------------------------------------------------------------

function setup_axes(Ax)
set(Ax,'tickdir','out','ticklength',[0.02 0.02],'linewidth',1,'fontsize',7,'xtick',0:20:120,'xlim',[0 100],'ylim',[0 1])
%set(gca,'ytick',0:0.1:1)
xlabel(Ax,'Distance along coastline','fontsize',10)

function setup_legend(ph)
legend(ph,'Gonochore','Fixed','Relative Size','Size Frequency','Location','northeast')

function plot_patches
patch([46, 46, 55, 55],[0, 1.0 1.0, 0],[-1, -1, -1, -1],[0.7 0.7 0.7],'edgecolor','none') ;

