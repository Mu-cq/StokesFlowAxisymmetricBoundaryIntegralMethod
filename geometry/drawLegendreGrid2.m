%draw a line having (x1,y1) as starting point and (x2,y2) as ending point
%and composed by N element

function [X,Y] = drawLegendreGrid2(theta,ampli)

    %second legendre polznomial
    P2 = legendre(2,cos(theta));
    P2 = P2(1,:);

    %perturbed radius
    R = (1+ampli*P2)';
    %R = alpha.*(1+ampli*P2);

    %build the circle (half)
    X = R.*cos(theta);
    Y = R.*sin(theta);

    %rescale in order to have volume equal to 4/3*pi
    V = axis_int_gauss_vect(X',Y');
    alpha = nthroot(4/3*pi/V,3);
    R = alpha*(1+ampli*P2)';

    %build the new circle (half)
    X = R.*cos(theta);
    Y = R.*sin(theta);

end