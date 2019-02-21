
function [X,Y,high,low]=draw_2modes(x1,y1,N,ampli1,ampli2)

    theta = 0:pi/N:pi;
    %second degree legendre polznomial
    P2 = legendre(2,cos(theta));
    P2 = P2(1,:);

    %third degree legendre polznomial
    P3 = legendre(3,cos(theta));
    P3 = P3(1,:);

    %perturbed radius
    R = (1+ampli1*P2+ampli2*P3);
    %R = alpha.*(1+ampli*P2);

    %build the circle (half)
    X = R.*cos(theta)+x1;
    Y = R.*sin(theta)+y1;

    %rescale in order to have volume equal to 4/3*pi
    V = axis_int_gauss(X,Y);
    alpha = nthroot(4/3*pi/V,3);
    R = alpha*(1+ampli1*P2+ampli2*P3);
    
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