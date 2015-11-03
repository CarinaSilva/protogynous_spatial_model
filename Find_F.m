function [F] = Find_F(FLEPtarg)

% Table lookup to find the desired value of F given a gonochore FLEP

load 'FLEPmat.mat' ;

Ftmp = abs(FLEPmat(:,2) - FLEPtarg) ;

index = find(Ftmp == min(Ftmp),1) ;

F = FLEPmat(index,1) ;
