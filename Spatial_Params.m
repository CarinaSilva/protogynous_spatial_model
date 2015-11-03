function Spatial_Params(PP,disp)

% Set up parameters for simulating network-persistent reserve scenarios

% Input: PP = number of spatial cells
%        disp = dispersal distance

% Load in life history parameters
load lifehistory_params.mat

T = 1000 ; % Number of time steps

N = nan(Amax,T,PP) ;
M = nan(Amax,T,PP) ;
SC = nan(Amax,T,PP) ;
fi = nan(Amax,T,PP) ;
Fi = nan(Amax,T,PP) ;

Settlers = nan(T,PP) ;
FertEggs = nan(PP,1);

% Dispersal matrix
Z = 1:PP ; % # cells
Z = repmat(Z,[PP,1]) ;  
mu = 0 ; % Mean

KK = normcdf(Z-Z'+0.5,mu,disp)-normcdf(Z-Z'-0.5,mu,disp); % Z-Z' finds distance, +/- 0.5 finds center of cell
KK2 = normcdf(Z-Z'+0.5+PP,mu,disp)-normcdf(Z-Z'-0.5+PP,mu,disp); % Adds in looping of cells
KK3 = normcdf(Z-Z'+0.5-PP,mu,disp)-normcdf(Z-Z'-0.5-PP,mu,disp); % Adds in looping other direction
KK4 = normcdf(Z-Z'+0.5+PP*2,mu,disp)-normcdf(Z-Z'-0.5+PP*2,mu,disp); % Adds in looping of cells
KK5 = normcdf(Z-Z'+0.5-PP*2,mu,disp)-normcdf(Z-Z'-0.5-PP*2,mu,disp); % Adds in looping other direction
Ktotal = KK+KK2+KK3+KK4+KK5 ; % Dji dispersal matrix
Ktotal = Ktotal/max(eig(Ktotal));

save spatial_params.mat


