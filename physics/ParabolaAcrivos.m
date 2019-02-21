% get ellipse shapes from linear theory Acrivos 1973

function [x,y] = ParabolaAcrivos(t,Ca)

x = 20*Ca^2*(1-2*t);
y = 1/4/Ca*(1-(x/20/Ca^2).^2);