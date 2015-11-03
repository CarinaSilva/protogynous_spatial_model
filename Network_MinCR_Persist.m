function Network_MinCR_Persist(Lf)

% This function finds the minimum proportion of coastline in reserves
% (C-sub-R in the manuscript) that allows population persistence for
% different reproductive life histories.

% Input is Lf, the mean size of entry to the fishery (range = 20 to 30)

FLEP = [0, 0.2]; % Fraction of lifetime egg production
S = {'GON','SC1','SC2','SC3'} ; % Life history scenarios GON
PHI = 1:1:20 ; % Male importance parameter
RS_init = 0.0 ; % initial guess at RS parameter (RS = C_R in the main text)

%  For network persistence, RW = 10, dist = 100
RW_master = 10;
disp = 100;

% run parameter creation file based on Lf
LifeHistory_Params(Lf); 
Gonochore_F_FLEP(Lf);

savename = strcat('spatialnetwork_minRS_Lf',num2str(Lf),'_persist','.mat');

MinRS = nan(length(S),length(PHI),length(FLEP));

for i = 1:length(S) 
    S{i}
for j = 1:length(PHI)
    PHI(j)
for f = 1:length(FLEP)
    
  
% find min fraction of coastline for persistence
RS= RS_init ;
Persist = false ;
Min = RS_init;
Max = 1;
while Max-Min>0.02
    
    if RS == 0
        RW = 0;
        PP = 10;
    else 
        RW = RW_master;
        PP = round(RW./RS); % number of patches
    end
    
% Run parameter creation file    
Spatial_Params(PP,disp) ; 
load('spatial_params.mat')

    % Get the appropriate value of FLEP
    F = Find_F(FLEP(f));
    
if FLEP(f) == 0;
    [~, ~, ~, Persist] = Spatial_Model(S{i},F,PHI(j),RW) ; % no repro outside MPAs
else
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, Persist] = Spatial_Model(S{i},F,PHI(j),RW) ;
end

    if Persist ;
        
        Max = RS; % this is the new max value
        Int = (Max-Min)/2;
        if Int < 0.1
            RS = RS - 0.01;
        else
            RS = RS - Int;
        end   
        
    else 
        Min = RS; % this is the new min value
        Int = (Max-Min)/2;
        if Int < 0.1
            RS = RS + 0.01;
        else
            RS = RS + Int;
        end
        
    end
    
    if Max-Min<=0.02
    MinRS(i,j,f) = Max; 
    end
 
    if RS > 1.0 % stopper in case we get crazy values
        MinRS(i,j,f) = NaN;
        Persist = true;
        Max = Min;
    end   

end % end while

end
end
end

save(savename) 