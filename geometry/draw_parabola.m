%draw a line having (x1,y1) as starting point and (x2,y2) as ending point
%and composed by N element

function [X,Y]=draw_parabola(x1,y1,x2,y2,N)

    %parabola coefficient
    a = (x2-x1)/(y2-y1)^2;
    
    %coordinates
    Y = y1:(y2-y1)/N:y2;
    X = a*(Y-y1).^2 + x1;
    
%     figure
%     plot(X,Y,'o-')
%     axis equal

end