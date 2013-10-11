
%plane-parallel with perfect rear mirror
function[af,ar] = PlaneParallel(alphaL, ~)
af = 1 - exp(-2*alphaL);
ar = 0;
end