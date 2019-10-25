function LifeHistory_Params(Lf)

% Specify life history parameters to be used in the model 

% Lf= comprimento (L) 

Amax = 20 ; % Number of age classes (max age is Amax-1 bc start counting at age 0)   % são 20 classes porque começa na idade 0 e termina em 20-1
            % No caso do Goraz vai da classe (ver atualizações no documento word com parâmetros a identificar no código matlab


a = 1 ; % constant in betacdf function (positive concavity)  % no paper diz que a concavidade é negativa

% Reproduction 
c = 7.04 ;  % Constant in fecundity relationship    % v no artigo tabela 2
e = 2.95 ;  % Exponent in allometric relationship   % tabela 2 
K = 0.000003 ;  % Slope of fert. fxn parameter     % 
X = 0.09 ;  % Intercept of fertilization fxn parameter %

% vonBert Growth
D = zeros(2,Amax) ; 
D(1,:) = (0:(Amax-1)) ;
D(2,1) = 8 ;
k = 0.05 ; %original shape
Linf = 90 ;
T0 = -1.875 ; % Age at size 0 (gives it 8cm at size 0)
D(2,:) = Linf.*(1-exp(-k.*((0:(Amax-1))-T0))) ;

% Recruitment 
%uA = 0.42
%uA = 0.28
uA = 0.35 ;  % Adult mortality
Surv = exp(-uA.*D(1,:)');
EggProd = c.*(D(2,:)').^e ;
% Stand-in maturity function to calculate LEP_unfished
Lm = 20 ; % Length at which 50% of fish mature 
q = 1 ;  % Shape parameter in maturity function 
% Probability of Maturity 
M_tmp = 1./(1+exp(-q.*((D(2,:)')-Lm)));  
% Stock-recruitment function
CRT = 0.25; % Critical value of LEP for persistence
LEP_unfished = ((sum(Surv.*EggProd.*M_tmp))/2);
Alpha = 1/(CRT*LEP_unfished);  % BevHolt params.
Beta = 1  ;

% Adult Survival
r =1;  % Steepness of selectivity curve 
isfished=nan(Amax,1) ; % Fishing mortality by age
isfished(:,1)=1./(1+exp(-r.*((D(2,:)')-Lf))) ;

save lifehistory_params.mat
clear 

% Parameters specific to each sex-change mode:
%--------------------------------------------------------------------------
% SC1 - Absolute Length
Lm = 20 ; % Length at which 50% of fish mature 
q = 1 ;  % Shape parameter in maturity function 
Lc = 30 ;  % Length at which 50% of fish change sex
p = 1 ; % Shape param. in sex change function

save SC1_params.mat
clear 

%--------------------------------------------------------------------------
% SC2 - Mean Length
Lm = 4 ; % Difference from the mean size at which prob. of maturity is 0.5
q = 1 ;  % Shape parameter in maturity function
Lc = 14 ; % Difference from the mean size at which prob. sex change is 0.5
p = 1 ;  % Shape parameter in sex change fxn

save SC2_params.mat
clear 

%--------------------------------------------------------------------------
% SC3 - Frequency of Smaller Individuals
Fm = 0.60 ;
%Fm = 0.50 ; % Frequency of smaller individuals where prob. of maturity is 0.5
q = 50 ;  % Shape parameter in maturity function
Fc = 0.67 ;
%Fc = 0.67 ; %********* %Fc = 0.67 ;  % Frequency of smaller mature indiv where prob of sex change is 0.5
p = 50 ;  % Shape parameter in sex change fxn

save SC3_params.mat


