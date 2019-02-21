%compute first and second derivative of a function computed with splines

function [fp,fpp] = derSplinesUnk(t,b,c,d)

fp = b+2*c.*t+3*d.*t.^2;
fpp = 2*c+6*d.*t;