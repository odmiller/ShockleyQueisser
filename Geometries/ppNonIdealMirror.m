
% T = transmission through rear mirror - assumed independent of angle
% ar = absorption for photons incident from the rear (i.e. at thermal
%   equilibrium)
% consistent with SQ, we take the front to have perfect ARC
function[af,ar] = ppNonIdealMirror(alphaL, ~, T, rearN)
af = T*(1-exp(-alphaL)) + (1-T)*(1-exp(-2*alphaL));

thetaC = asin(sqrt(1/rearN^2));
ar = @(theta) rearN^2*T*(1-exp(-alphaL/cos(theta)))*(theta<=thetaC) ...
    + rearN^2*T*(1-exp(-2*alphaL/cos(theta)))./(1-(1-T)*exp(-2*alphaL/cos(theta)))*(theta>thetaC) ;

end