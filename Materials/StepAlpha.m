
% step-function alpha
% analytic approximation to step-function, can make kStep arbitrarily large
function[alphaL, n] = StepAlpha(E, Eg, alphaL0)
    kStep = 150; % roughly that of GaAs
    alphaL = alphaL0 ./( 1 + exp(-2*kStep*(E-Eg)));
    n = 3.5;
end