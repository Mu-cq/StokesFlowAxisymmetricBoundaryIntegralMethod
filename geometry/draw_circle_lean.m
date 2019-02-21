
%draw a line having (x1,y1) as starting point and (x2,y2) as ending point
%and composed by N element

function [X,Y,high,low]=draw_circle_lean(x1,y1,N,DELTA)

    %indipendent variable
    theta = 0:pi/N:pi;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %semi axis from D and volume of sphere with radius=1
    b = ((1.0-DELTA)/(1.0+DELTA))^(1/3);
    a = 1.0/b^2;


    %build drop, volume is the one of a sphere of radius 1
    X = a*cos(theta) + x1;
    Y = b*sin(theta) + y1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ds = sqrt((X(1:end-1)-X(2:end)).^2+(Y(1:end-1)-Y(2:end)).^2);
    high = max(ds);
    low = min(ds);

end