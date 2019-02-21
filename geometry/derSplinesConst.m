%compute first and second derivative of a function computed with splines

function [fp,fpp] = derSplinesConst(b,c,d)

fp = @(t) b+2*c*t+3*d*t^2;
fpp = @(t) 2*c+6*d*t;

fp = fp(0.5);
fpp = fpp(0.5);