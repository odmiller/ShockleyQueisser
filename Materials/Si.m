
% crystalline Silicon
% alpha = absorption coefficient, n = refractive index
% note units of alphaBeta
function[alphaL, n, alphaBeta] = Si(~, L)
    load('alphaSi.mat');
    alphaL = L * alphaNew;
    n = 3.5; % 4n^2
    
    C=3.88E-31;  % Auger coefficient
    ni=1.45E10;
    alphaBeta = [C*ni^3*L 3/2]; % Auger
end