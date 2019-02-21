%compute cubic NATURAL spline coefficient

function [ax, bx, cx, dx, ay, by, cy, dy] = my_spline_periodic (x, y)
    
    A = diag(4*ones(1,numel(x)));
    A(1,end-1) = 1;
    A(end,2) = 1;
    %A(1,1) = 2;
    %A(end,end) = 2;
    A = A+diag(1*ones(1,numel(x)-1),1);
    A = A+diag(1*ones(1,numel(x)-1),-1);
    
    deltaX = x(3:end)-x(1:end-2);
    %deltaX = [x(2)-x(end) deltaX x(1)-x(end-1)];
    deltaX = [x(2)-x(end-1) deltaX x(2)-x(end-1)];
    
    deltaY = y(3:end)-y(1:end-2);
    %deltaY = [y(2)-y(end) deltaY y(1)-y(end-1)];
    deltaY = [y(2)-y(end-1) deltaY y(2)-y(end-1)];
    
    Dx = A\(3*deltaX');
    Dy = A\(3*deltaY');
    
    ax = x(1:end-1);
    bx = Dx(1:end-1)';
    cx = 3*(x(2:end)-x(1:end-1))-2*Dx(1:end-1)'-Dx(2:end)';
    dx = 2*(x(1:end-1)-x(2:end))+Dx(1:end-1)'+Dx(2:end)';
    
    ay = y(1:end-1);
    by = Dy(1:end-1)';
    cy = 3*(y(2:end)-y(1:end-1))-2*Dy(1:end-1)'-Dy(2:end)';
    dy = 2*(y(1:end-1)-y(2:end))+Dy(1:end-1)'+Dy(2:end)';
    
%     figure(10)
%     plot(x,y,'o')
%     hold on
%     t = 0:1/1000:1;
%     for i = 1:numel(x)-1
%         xx = ax(i)+bx(i)*t+cx(i)*t.^2+dx(i)*t.^3;
%         yy = ay(i)+by(i)*t+cy(i)*t.^2+dy(i)*t.^3;
%         plot(xx,yy,'r-')
%     end
%     axis equal
%     hold off
    
    %check continuity of first derivative

return