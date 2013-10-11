

% GaAs
% alpha = absorption coefficient, n = refractive index
% L = thickness, in cm
% model as used in 
%	Miller et. al., IEEE Jour. of Photovolt. 2, 303 (2012)
function[alphaL, n, alphaBeta] = GaAs(E, L) % GaAs parameters
    Eg = 1.4255;
    E0 = 0.0067;
    alpha0 = 8000;
    Eprime = 0.140; % best fit to data (src)
    deltaE = E - Eg;
    alphaL = alpha0 * L ...
        *( exp(deltaE/E0) .* (E<=Eg) + (1 + deltaE/Eprime) .* (E>Eg) );
    n = 3.5;
    
    C=7E-30;    % Auger
    ni=1.79E6;
    alphaBeta = [C*ni^3 3/2]; 
end
