%find Parabola coeefiicient passing trough two point and symmetric to the
%axis

function [a,b,c] = ParabolaCoeff(x1,y1,x2,y2)

    a = (y1-y2)./(x1.^2-x2.^2);
    b = 0;
    c = y2-a*x2.^2;

end