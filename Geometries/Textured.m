
% textured with perfect rear mirror
function[af,ar] = Textured(alphaL, n)
af = alphaL ./ (alphaL + 1/(4*n^2));
ar = 0;
end