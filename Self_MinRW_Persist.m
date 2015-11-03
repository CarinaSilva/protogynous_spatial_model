function Self_MinRW_Persist(Lf)

% This function finds the minimum proportion of coastline in reserves
% (C-sub-R in the manuscript) that allows population persistence for
% different reproductive life histories.

% Input is Lf, the mean size of entry to the fishery (range = 20 to 30)

FLEP = [0, 0.2]; % Fraction of lifetime egg production
S = {'GON','SC1','SC2','SC3'} ; % Life history scenarios GON
PHI = 1:1:20 ; % Male importance parameter
RW_init = 0.0 ; % initial guess at RW parameter (RW = R_W in the main text)

% For self persistence, RS = 0.1, dist = 10;
RS_master = 0.1;
disp = 10;

% run parameter creation file based on Lf
LifeHistory_Params(Lf); 
Gonochore_F_FLEP(Lf);

savename = strcat('spatialself_minRW_Lf',num2str(Lf),'_persist','.mat'); 

MinRW = nan(length(S),length(PHI),length(FLEP));

for i = 1:length(S)
    S{i}
for j = 1:length(PHI)
    PHI(j)
    
    if strcmp(S{i},'GON') && j > 1
        MinRW(i,j,:) = MinRW(i,1,:);
    else
    
for f = 1:length(FLEP)
    
% find min reserve width for persistence
RW = RW_init ;
Persist = false ;
Min = RW_init;
Max = 100;
while Max-Min > 1
    
    if RW == 0
        RS = 0;
        PP = 10;
    else
        RS = RS_master;
        PP = round(RW/RS); % number of patches
    end
    
Spatial_Params(PP,disp) ; % Run parameter creation file
load('spatial_params.mat')

    % Get the new value of FLEP
    F = Find_F(FLEP(f));
    
if FLEP(f) == 0;
    [~, ~, ~, Persist] = Spatial_Model(S{i},F,PHI(j),RW) ; % no repro outside MPAss
else
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, Persist] = Spatial_Model(S{i},F,PHI(j),RW) ;
end

    if Persist
        
        Max = RW;
        Int = round((Max-Min)/2);
        if Int < 10
            RW = RW - 1;
        else
            RW = RW - Int;
        end
        
    else
        Min = RW;
        Int = round((Max-Min)/2);
        if Int < 10;
            RW = RW + 1;
        else
            RW = RW + Int;
        end
    
    end
        if Max-Min<2
        MinRW(i,j,f) = Max ;
        end

if RW > 200 % stopper in case we get crazy values
MinRW(i,j,f) = NaN;
Persist = true;
Max=Min;
end

end % end while

end % loop over FLEP
end % end if GON & PHI(1)
end
end

save(savename)