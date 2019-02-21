%compute the coefficient of the equation for a circoference given three
%points

function [R, a, b, c] = CoeffCirc(x1, y1, x2, y2, x3, y3)

    b = ((x2*x2+y2*y2-x3*x3-y3*y3)*(x2-x1)-(x1*x1+y1*y1-x2*x2-y2*y2)*(x3-x2))/((y1-y2)*(x3-x2)-(y2-y3)*(x2-x1));
    a = (x2*x2+y2*y2-x3*x3-y3*y3+b*(y2-y3))/(x3-x2);
    c = -(x3*x3+y3*y3+a*x3+b*y3);
    
    R = sqrt(-c+a/2*a/2+b/2*b/2);

end