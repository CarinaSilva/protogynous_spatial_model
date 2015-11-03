function Spatial_Struct(PHI,FLEP,Type)

% This function organizes a set of model runs for a range of Phi and FLEP
% values. This is mostly used for plotting up spatial results for different
% scenarios.

% Inputs: PHI (male importance Beta parameter)
%         FLEP (fraction of lifetime egg production in gonochore -
%         specifies equivalent F value
%         Type: type of persistence possible: network or self (string)

% Scenario for figures:
switch Type
    case 'network'
disp = 100; % dispersal distance
RW = 10; % Reserves are too small for self-persistence (dispersal distance = 100)
RS = 0.25; % Reserve fraction = CRT, so network persistence is possible.
PP = round(RW./RS); % number of patches

    case 'self'
disp = 10; % dispersal distance
RW = 10; % Reserves are just large enough for self-persistence (dispersal distance = 10)
RS = 0.1; % Reserve fraction < CRT, so network persistence is not possible.
PP = round(RW./RS); % number of patches

end % end switch over Type

Lf = 20; % Size of entry to the fishery
LifeHistory_Params(Lf);
Gonochore_F_FLEP(Lf);

Spatial = struct([]); % Create empty structure
Spatial_Params(PP,disp) ; % Run parameter creation file
load('spatial_params.mat')

datename = 'Oct2015';

savename = strcat('spatial',Type,'_',datename,'_PHI',num2str(PHI),'_FLEP',num2str(FLEP),'_baseline','.mat') ;
 
S = {'GON','SC1','SC2','SC3'}; % Scenarios

% Convert FLEP to F
Gonochore_F_FLEP(Lf)
F = Find_F(FLEP);

% Loop over simulations
for s = 1:length(S) ;
    
for f = 1:length(F)
    Spatial(1).(S{s}).F(f).F=F(f) ; % store metadata

for phi = 1:length(PHI)
    Spatial(1).(S{s}).F(f).PHI(phi).PHI=PHI(phi) ; % store metadata
  
for rs = 1:length(RS)
    
    % store metadata
    Spatial(1).(S{s}).F(f).PHI(phi).RS(rs).RS=RS(rs) ;
    Spatial(1).(S{s}).F(f).PHI(phi).RS(rs).RW=RW(rs) ;

    [AgeDist, Ntotal, Biomass, Persist, FBiom, MBiom, MBiomPatch, Sexratio, SexratioNum, ProbSC, ProbMat, Persist2,BiomPatch] = Spatial_Model(S{s},F(f),PHI(phi),RW(rs)) ;
    
    % Save results:
    Spatial(1).(S{s}).F(f).PHI(phi).RS(rs).AgeDist = AgeDist ;
    Spatial(1).(S{s}).F(f).PHI(phi).Ntotal(rs) = Ntotal ;
    %Spatial(1).(S{s}).F(f).PHI(phi).Ytotal(rs) = Ytotal ; 
    %Spatial(1).(S{s}).F(f).PHI(phi).RS(rs).Sratio = Sratio ;
    Spatial(1).(S{s}).F(f).PHI(phi).Biomass(rs) = Biomass ; 
    Spatial(1).(S{s}).F(f).PHI(phi).Persist(rs) = Persist;
    Spatial(1).(S{s}).F(f).PHI(phi).FBiom(rs) = FBiom; %female biomass
    Spatial(1).(S{s}).F(f).PHI(phi).MBiom(rs) = MBiom; %male biomass
    Spatial(1).(S{s}).F(f).PHI(phi).RS(rs).MBiomPatch = MBiomPatch ; %male biomass per patch
    Spatial(1).(S{s}).F(f).PHI(phi).RS(rs).Sexratio = Sexratio ; %biomass sex ratio
    Spatial(1).(S{s}).F(f).PHI(phi).RS(rs).SexratioNum = SexratioNum ; %numerical sex ratio
    Spatial(1).(S{s}).F(f).PHI(phi).RS(rs).ProbSC = ProbSC  ; %probability of sex change
    Spatial(1).(S{s}).F(f).PHI(phi).RS(rs).ProbMat = ProbMat  ; %probability of maturation
    Spatial(1).(S{s}).F(f).PHI(phi).Persist2(rs) = Persist2; % persistence, accounting for reproduction elsewhere
    Spatial(1).(S{s}).F(f).PHI(phi).RS(rs).BiomPatch = BiomPatch; % total biomass in each patch
end
end
end
end

save(savename)



