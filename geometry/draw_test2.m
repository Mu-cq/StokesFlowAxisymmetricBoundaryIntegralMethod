%draw test inner

function [x,y] = draw_test2(x1,y1,L,R,N)
    
    %create the empty matrix
    theta1 = 0:pi/(N/4)/2:pi/2;
    theta2 = pi/2:pi/(N/4)/2:pi;

    %build the semicircles (half)
    X1 = R.*cos(theta1)+x1+L/2;
    Y1 = R.*sin(theta1)+y1;
    X3 = R.*cos(theta2)+x1-L/2;
    Y3 = R.*sin(theta2)+y1;
    
    %build straigth part
    [X2,Y2]=drawline2(X1(end),Y1(end),X3(1),Y3(1),N/2);
    %[X4,Y4]=drawline2(X3(end),Y3(end),X1(1),Y1(1),N/4);
    
    x = [X1(1:end-1) X2(1:end-1) X3];
    y = [Y1(1:end-1) Y2(1:end-1) Y3];
    
%     figure
%     plot(x,y,'o-')
%     axis equal
%     xlabel('x')
%     ylabel('y')
    

end