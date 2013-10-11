
% crystalline Silicon
% alpha = absorption coefficient, n = refractive index
% note units of alphaBeta
% Si data digitized (not great) from 
%	Tiedje et. al., IEEE Trans. on Electron Devices 31, 711 (1984)
% Auger taken from same paper
% A more correct model of the Auger coefficient can be found in 
%	Kerr et. al., PVSC 2002, "Lifetime and efficiency limits of crystalline
%	silicon solar cells"
function[alphaL, n, alphaBeta] = Si(~, L)
    load('alphaSi.mat');
    alphaL = L * alphaNew;
    n = 3.5; % 4n^2
    
    C=3.88E-31;  % Auger coefficient
    ni=1.45E10;
    alphaBeta = [C*ni^3 3/2]; % Auger
end