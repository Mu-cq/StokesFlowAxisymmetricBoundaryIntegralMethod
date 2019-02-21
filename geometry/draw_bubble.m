%draw a line having (x1,y1) as starting point and (x2,y2) as ending point
%and composed by N element

function [X,Y,high,low,alpha]=draw_bubble(x1,y1,N,DELTA,R)

    theta = 0:pi/N:pi;

    a = 1;
    b = (1-DELTA)/(1+DELTA);

    %build the circle (half)
    X = a.*cos(theta)+x1;
    Y = b.*sin(theta)+y1;

    %rescale in order to have volume equal to 4/3*pi*R^3
    V = axis_int_gauss(X,Y);
    alpha = nthroot(4/3*pi*R^3/V,3);
    a = alpha*a;
    b = alpha*b;

    %build the new circle (half)
    X = a.*cos(theta)+x1;
    Y = b.*sin(theta)+y1;

    ds = sqrt((X(1:end-1)-X(2:end)).^2+(Y(1:end-1)-Y(2:end)).^2);
    high=max(ds);
    low=min(ds);

end