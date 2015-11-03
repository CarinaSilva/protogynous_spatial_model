function [AgeDist, Ntotal, Biomass, Persist, FBiom, MBiom, MBiomPatch, Sexratio, SexratioNum, ProbSC, ProbMat, Persist2,BiomPatch] = Spatial_Model(S,F,PHI,RW)

load spatial_params.mat
load lifehistory_params.mat

% Set up spatial domain
XX = [ones(1,PP-RW),zeros(1,RW)];
patchfished = XX./( (PP-RW)/PP);% Which patches are fished 


% Initialize population          
for pp = 1:PP

% Create Leslie matrix
     Z = uA ; 
     surv = exp(-Z);   % Instantaneous mortality rate... e^-M
     Surv = repmat(surv,[Amax,1]) ;
    
     N(:,1,pp) = [1; cumprod(Surv(1:end-1))]; % Put in one recruit & let it survive through age classes.
      
switch S 

%--------------------------------------------------------------------
    case 'GON'  %Gonochores
    
%-------------------------------------------------------------------------    
case 'SC1'   % Absolute Length
load SC1_params.mat

% Probability of maturity
M(:,1,pp) = 1./(1+exp(-q.*((D(2,:)')-Lm)));  
% Probability of sex Change
SC(:,1,pp) = 1./(1+exp(-p.*((D(2,:)')-Lc))) ; % Prob. of being male
TB(1,pp) = sum((N(:,1,pp).*M(:,1,pp).*D(2,:)'.^3)) ; % Total mature biomass
MB(1,pp) = sum((N(:,1,pp).*SC(:,1,pp).*M(:,1,pp).*D(2,1:end)'.^3)) ; % Total mature male biomass

%------------------------------------------------------------------------
case 'SC2'   % Mean Length
load SC2_params.mat

% Probability of maturity
M(:,1,pp) = 1./(1+exp(-q.*((D(2,:)'-((sum(D(2,:)'.*N(:,1,pp))/(sum(N(:,1,pp))))+Lm))))) ;
% Probability of sex Change
SC(:,1,pp) = 1./(1+exp(-p.*((D(2,:)'-((sum(D(2,:)'.*N(:,1,pp))/(sum(N(:,1,pp))))+Lc))))) ; %prob of being male
TB(1,pp) = sum((N(:,1,pp).*M(:,1,pp).*D(2,:)'.^3)) ; % Total mature biomass
MB(1,pp) = sum((N(:,1,pp).*SC(:,1,pp).*M(:,1,pp).*D(2,1:end)'.^3)) ; % Total mature male biomass

%-------------------------------------------------------------------------
case 'SC3'   % Frequency of Smaller Individuals
load SC3_params.mat

% Probability of Maturity
fi(:,1,pp) = cumsum(N(:,1,pp)./sum(N(:,1,pp)));
fi(:,1,pp) = [0; fi(1:end-1,1,pp)];
M(:,1,pp) = 1./(1+exp(-q.*(fi(:,1,pp)-Fm))) ;

% Probability of Sex Change
Fi(:,1,pp) = cumsum((N(:,1,pp).*M(:,1,pp))./sum((N(:,1,pp).*M(:,1,pp)))) ;
Fi(:,1,pp) = [0; Fi(1:end-1,1,pp)] ;
SC(:,1,pp) = 1./(1+exp(-p.*(Fi(:,1,pp)-Fc)))  ; 
TB(1,pp) = sum((N(:,1,pp).*M(:,1,pp).*D(2,:)'.^3)) ; % Total mature biomass
MB(1,pp) = sum((N(:,1,pp).*SC(:,1,pp).*M(:,1,pp).*D(2,1:end)'.^3)) ; % Total mature male biomass

end % end switch
end
 
% Start running model w/ initial population params.   
for pp = 1:PP    
    
     Ftmp(pp) = patchfished(pp).*F ;
     Z = (uA + Ftmp(pp).*isfished) ;
     surv = exp(-Z);   % Instantaneous mortality rate... e^-M
     N(:,1,pp) = cumprod([1; surv(1:end-1)]) ;

     Atmp = diag(surv(1:end-1)) ; 
     Atmp = [Atmp,zeros(size(Atmp,1),1)] ; 
     Atmp = [zeros(1,size(Atmp,2));Atmp] ;
     A(:,:,pp) = Atmp ;    
    
end

for t = 2:T
    for pp = 1:PP
       N(:,t,pp) = A(:,:,pp)*N(:,t-1,pp) ;
   
%--------------------------------------------------------------------------
% Reproduction switch
switch S

case 'GON'   % Gonochores
        
% Reproduction        
E =(EggProd.*M_tmp.*N(:,t-1,pp)./2) ;
FertEggs(pp)=sum(E);
B(t,pp) = 0.5 ;
Pf(t,pp) = 1.0; %fillers for sex changer calculations
TB(t,pp) = 1.0;
MB(t,pp) = 1.0;
B2(t,pp) = 0.5 ;

case {'SC1','SC2','SC3'} ;  % Sex-changing
    
% Reproduction
E(t,pp) = sum((c.*D(2,:)'.^e).*N(:,t-1,pp).*(1-SC(:,t-1,pp)).*M(:,t-1,pp)) ; % Egg production
TB(t,pp) = sum((N(:,t-1,pp).*M(:,t-1,pp).*D(2,:)'.^3)) ; % Total mature biomass
MB(t,pp) = sum((N(:,t-1,pp).*SC(:,t-1,pp).*M(:,t-1,pp).*D(2,1:end)'.^3)) ; % Total mature male biomass
B(t,pp) = MB(t-1,pp)/TB(t-1,pp) ;  % Proportion of mature male biomass (Biomass sex ratio)
B2(t,pp)= sum(N(:,t-1,pp).*M(:,t-1,pp).*SC(:,t-1,pp))./sum(N(:,t-1,pp).*M(:,t-1,pp)); % Numerical sex ratio 
Pf(t,pp) = (betacdf(B(t,pp),a,PHI)) ; % Proportion of eggs fertilized
FertEggs(pp) = E(t,pp)*Pf(t,pp) ;

end % end reproduction switch
    end 
    
for pp = 1:PP  
    
% Recruits entering A0
Settlers(t,:)=Ktotal*FertEggs(:) ;
N(1,t,pp) = ((Alpha.*Settlers(t,pp))./(1+Alpha./Beta.*Settlers(t,pp)))  ; % BevHolt

switch S 

%--------------------------------------------------------------------
case 'GON'  %Gonochores
%-------------------------------------------------------------------------    
case 'SC1'   % Absolute Length
load SC1_params.mat
% Probability of maturity
M(:,t,pp) = 1./(1+exp(-q.*((D(2,:)')-Lm)));  
% Probability of sex Change
SC(:,t,pp) = 1./(1+exp(-p.*((D(2,:)')-Lc))) ; % Prob. of being male

%------------------------------------------------------------------------
case 'SC2'   % Mean Length
load SC2_params.mat
% Probability of maturity
M(:,t,pp) = 1./(1+exp(-q.*((D(2,:)'-((sum(D(2,:)'.*N(:,t,pp))/(sum(N(:,t,pp))))+Lm))))) ;
% Probability of sex Change
SC(:,t,pp) = 1./(1+exp(-p.*((D(2,:)'-((sum(D(2,:)'.*N(:,t,pp))/(sum(N(:,t,pp))))+Lc))))) ; %prob of being male

%-------------------------------------------------------------------------
case 'SC3'   % Frequency of Smaller Individuals
load SC3_params.mat
% Probability of Maturity
fi(:,t,pp) = cumsum(N(:,t,pp)./sum(N(:,t,pp)));
fi(:,t,pp) = [0; fi(1:end-1,t,pp)];
M(:,t,pp) = 1./(1+exp(-q.*(fi(:,t,pp)-Fm))) ;
% Probability of Sex Change
Fi(:,t,pp) = cumsum((N(:,t,pp).*M(:,t,pp))./sum((N(:,t,pp).*M(:,t,pp)))) ;
Fi(:,t,pp) = [0; Fi(1:end-1,t,pp)] ;
SC(:,t,pp) = 1./(1+exp(-p.*(Fi(:,t,pp)-Fc)))  ; 

end % end switch
% To calculate yield:
%Y(t,pp) = sum((N(1:end-1,t-1,pp)-N(2:end,t,pp)).*(Ftmp(pp).*isfished(1:end-1))./(uA+(Ftmp(pp).*isfished(1:end-1))).*(D(2,1:end-1)').^3) ;
  
end
end

%CHECK FOR PERSISTENCE - SCORCHED EARTH
for pp = 1:PP

        Ftmp(pp) = patchfished(pp).*F ;
        Z = (uA + Ftmp(pp).*isfished) ;
        surv = exp(-Z);   % Instantaneous mortality rate... e^-M
        NN(:,pp) = cumprod([1; surv(1:end-1)]) ;
        
        if patchfished(pp)
            NN(:,pp) = 0;
        end
        
        switch S
            case 'GON'
                EE(pp) =sum(EggProd.*M_tmp.*NN(:,pp)./2) ;
                FE(pp)=(EE(pp));

            case {'SC1','SC2','SC3'}
                EE(pp) = sum((c.*D(2,:)'.^e).*NN(:,pp).*(1-SC(:,end,pp)).*M(:,end,pp));
                 FE(pp) = EE(pp).*Pf(end,pp);
       end
       
end

%CHECK FOR PERSISTENCE - REPRO OUTSIDE RESERVE
for pp = 1:PP

        Ftmp(pp) = patchfished(pp).*F ;
        Z = (uA + Ftmp(pp).*isfished) ;
        surv = exp(-Z);   % Instantaneous mortality rate... e^-M
        NNN(:,pp) = cumprod([1; surv(1:end-1)]) ;
        
        switch S
            case 'GON'
                EE2(pp) =sum(EggProd.*M_tmp.*NNN(:,pp)./2) ;
                FE2(pp)=(EE2(pp));

            case {'SC1','SC2','SC3'}
                EE2(pp) = sum((c.*D(2,:)'.^e).*NNN(:,pp).*(1-SC(:,end,pp)).*M(:,end,pp));
                 FE2(pp) = EE2(pp).*Pf(end,pp);
       end
       
end


%----------------------------------------------------------------------
% Outputs
AgeDist = N(:,end,:) ;
Ntotal = sum(sum(N(:,end,:))) ; % Ntotal for each F value
%Ytotal = sum(Y(end,:)) ;
Biomassatlength = D(2,:).^3;
Popsize = squeeze(N(:,end,:))' ;

Biomass = mean(Popsize*Biomassatlength') ;
BiomPatch = Popsize*Biomassatlength';
FBiom = (sum(TB(end,:))) - (sum(MB(end,:))) ;
MBiom = sum(MB(end,:)) ;
MBiomPatch = MB(end,:) ;

Sexratio = B(end,:) ; % Sex Ratio
SexratioNum = B2(end,:) ; % Numerical Sex ratio
ProbSC = SC(:,end,:) ;
ProbMat = M(:,end,:) ;

% Check for Persistence - SCORCHED EARTH
EPR = FE ; % eggs per recruit in the last timestep
C = repmat(EPR(:)',[PP,1]).*Alpha ; 
if any(isnan(C))
    Persist = 0;
else
[V,L] = eig(C.*Ktotal) ;
L = diag(L) ;
L1 = L(L==max(L))  ;%dom. eigenvalue
V1 = V(:,L==max(L)) ; %dom. eigenvector

if L1 > 1
   Persist = 1 ;
else Persist = 0 ;

end
end

%----------------------------------------------------------------------

% Check for Persistence - REPRO OUTSIDE RESERVE
EPR2 = FE2 ; % eggs per recruit in the last timestep
C = repmat(EPR2(:)',[PP,1]).*Alpha ; 
if any(isnan(C))
    Persist2 = 0;
else
[V,L] = eig(C.*Ktotal) ;
L = diag(L) ;
L1 = L(L==max(L))  ;%dom. eigenvalue
V1 = V(:,L==max(L)) ; %dom. eigenvector

if L1 > 1
   Persist2 = 1 ;
else Persist2 = 0 ;

end
end


