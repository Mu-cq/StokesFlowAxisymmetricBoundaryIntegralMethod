%draw test inner

function [x,y] = draw_test_choose2(x1,y1,R,L,N,M,theta)
    
    %for circles
    theta1 = pi/2:-pi/M:-pi/2 + theta;
    theta2 = 3*pi/2:-pi/M:pi/2 + theta;

    %build the semicircles (half)
    X1 = R.*cos(theta1) + x1 + L/2*cos(theta);
    Y1 = R.*sin(theta1) + y1 + L/2*sin(theta);
    X3 = R.*cos(theta2) + x1 - L/2*cos(theta);
    Y3 = R.*sin(theta2) + y1 - L/2*sin(theta);
    
    %build straigths parts
    [X2,Y2]=drawline2(X1(end),Y1(end),X3(1),Y3(1),N);
    [X4,Y4]=drawline2(X3(end),Y3(end),X1(1),Y1(1),N);
    
    x = [X1(1:end-1) X2(1:end-1) X3(1:end-1) X4];
    y = [Y1(1:end-1) Y2(1:end-1) Y3(1:end-1) Y4];
    
%     figure
%     plot(x,y,'o-')
%     axis equal
%     xlabel('x')
%     ylabel('y')
    

end