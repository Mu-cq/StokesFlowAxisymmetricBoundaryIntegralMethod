%draw test inner

function [x,y] = draw_test3(x1,y1,L,R1,R2,N,theta)
    
    %create the empty matrix
    theta1 = 0:(pi/2+theta)/(N/4):pi/2+theta;
    theta2 = pi/2+theta:(pi/2-theta)/(N/4):pi;

    %build the semicircles (half)
    X1 = R1.*cos(theta1)+x1+L/2;
    Y1 = R1.*sin(theta1)+y1;
    X3 = R2.*cos(theta2)+x1-L/2;
    Y3 = R2.*sin(theta2)+y1;
    
    %build straigth part
    [X2,Y2]=drawline2(X1(end),Y1(end),X3(1),Y3(1),N/2);
    
    x = [X1(1:end-1) X2(1:end-1) X3];
    y = [Y1(1:end-1) Y2(1:end-1) Y3];

end