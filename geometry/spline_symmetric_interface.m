%compute cubic spline coefficient

function [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric_interface(x, y)

    A = zeros(numel(x));
    
    A(2:end-1,2:end-1) = diag(4*ones(1,numel(x)-2));
    A(1,1) = 2;
    A(end,end) = 2;
    A = A+diag(1*ones(1,numel(x)-1),1);
    A = A+diag(1*ones(1,numel(x)-1),-1);
    
    deltaX = x(3:end)-x(1:end-2);
    deltaX = [x(2)-x(1) deltaX x(end)-x(end-1)];
    
    deltaY = y(3:end)-y(1:end-2);
    deltaY = [y(2)-y(1) deltaY y(end)-y(end-1)];
    
    Dx = A\(3*deltaX');
    
    %Symmetry conditions
    deltaY(1) = 0;
    deltaY(end) = 0;
    A(1,1) = 1;
    A(1,2) = 0;
    A(end,end-1) = 0;
    A(end,end) = 1;
    
    Dy = A\(3*deltaY');
    
    ax = x(1:end-1);
    bx = Dx(1:end-1)';
    cx = 3*(x(2:end)-x(1:end-1))-2*Dx(1:end-1)'-Dx(2:end)';
    dx = 2*(x(1:end-1)-x(2:end))+Dx(1:end-1)'+Dx(2:end)';
    
    ay = y(1:end-1);
    by = Dy(1:end-1)';
    cy = 3*(y(2:end)-y(1:end-1))-2*Dy(1:end-1)'-Dy(2:end)';
    dy = 2*(y(1:end-1)-y(2:end))+Dy(1:end-1)'+Dy(2:end)';

return