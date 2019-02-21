%compute first and second derivative of a function computed with splines

function [fp,fpp] = derSplines(b,c,d)

fp = [b b(end)+2*c(end)+3*d(end)];
fpp = [2*c 2*c(end)+6*d(end)];