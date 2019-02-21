%draw a line having (x1,y1) as starting point and (x2,y2) as ending point
%and composed by N element

function [X,Y,high,low,alpha,theta]=draw_circle_clean(x1,y1,N,ampli)

    theta = 0:pi/N:pi;
    %second legendre polznomial
    P2 = legendre(2,cos(theta));
    P2 = P2(1,:);

    %perturbed radius
    R = (1+ampli*P2);
    %R = alpha.*(1+ampli*P2);

    %build the circle (half)
    X = R.*cos(theta)+x1;
    Y = R.*sin(theta)+y1;

    %rescale in order to have volume equal to 4/3*pi
    V = axis_int(X,Y);
    alpha = nthroot(4/3*pi/V,3);
    R = alpha*(1+ampli*P2);

    %build the new circle (half)
    X = R.*cos(theta)+x1;
    Y = R.*sin(theta)+y1;

    ds = sqrt((X(1:end-1)-X(2:end)).^2+(Y(1:end-1)-Y(2:end)).^2);
    high=max(ds);
    low=min(ds);

%     figure
%     plot(X,Y,'o-',X,-Y,'o-b')
%     axis equal

end