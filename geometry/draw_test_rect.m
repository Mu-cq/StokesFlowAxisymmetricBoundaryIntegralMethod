%draw test inner

function [x,y] = draw_test_rect(x1,y1,R,L,N,M)
    
    %create the empty matrix
%     theta1 = -pi/2:pi/M:pi/2;
%     theta2 = pi/2:pi/M:3*pi/2;
% 
%     %build the semicircles (half)
%     X1 = R.*cos(theta1)+x1+L/2;
%     Y1 = R.*sin(theta1)+y1;
%     X3 = R.*cos(theta2)+x1-L/2;
%     Y3 = R.*sin(theta2)+y1;
    
    %build straigths parts
    [X1,Y1]=drawline2(L/2+x1,-R/2+y1,L/2+x1,R/2+y1,M);
    [X2,Y2]=drawline2(L/2+x1,R/2+y1,-L/2+x1,R/2+y1,N);
    [X3,Y3]=drawline2(-L/2+x1,R/2+y1,-L/2+x1,-R/2+y1,M);
    [X4,Y4]=drawline2(-L/2+x1,-R/2+y1,L/2+x1,-R/2+y1,N);
    
    x = [X1(1:end-1) X2(1:end-1) X3(1:end-1) X4];
    y = [Y1(1:end-1) Y2(1:end-1) Y3(1:end-1) Y4];
    
%     figure
%     plot(x,y,'o-')
%     axis equal
%     xlabel('x')
%     ylabel('y')
    

end