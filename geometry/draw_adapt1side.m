%draw a line having (x1,y1) as starting point and (x2,y2) as ending point
%and composed by N element

function [X,Y]=draw_adapt1side(x1,y1,x2,y2,N)

    %parameter space
    eta = x1:(x2-x1)/N:x2;
    beta = y1:(y2-y1)/N:y2;

    %physical space
    X = x1 + (x2-x1)*exp(-(eta-x1)/(x2-x1)+1) - 1;
    Y = y1 + (y2-y1)*exp(-(beta-y1)/(y2-y1)+1) - 1;

    figure
    plot(X,Y,'o-')
    axis equal
    grid on

end