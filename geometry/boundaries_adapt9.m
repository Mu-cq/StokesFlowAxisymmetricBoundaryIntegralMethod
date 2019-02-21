%create two vectors containing the boundaries nodes (first point and last point of the element) once you have given the radius R, the
%lenght L of the pipe, the number of element in each region (n for the
%outlet, m for the wall, j for the inlet).

function [a,b]=boundaries_adapt9(PARAM)

    m = PARAM.m;    q = PARAM.q;
    m1 = PARAM.m1;     m2 = PARAM.m2;    
    l = PARAM.L;
    theta = PARAM.theta;
    R = PARAM.R;
    thickness = PARAM.thickness;

    %build lower panel
    phi = 0:1/m1:1;
    xLow = (1-phi)*l;
    yLow = -thickness*ones(1,numel(phi)) + R;
    
    %build upper panel
    phi = 0:1/m1:1;
    xUp = phi*l;
    yUp = thickness*ones(1,numel(phi)) + R;
    
    %build left circle
    phi = 0:1/m2:1;
    s = pi/2+pi*phi;
    xLeft = thickness*cos(-s);
    yLeft = thickness*sin(-s) + R;
    
    %build left circle
    phi = 0:1/m2:1;
    s = 3*pi/2+pi*phi;
    xRight = thickness*cos(-s) + l;
    yRight = thickness*sin(-s) + R;
    
    %paste wall
    a(1:m+1) = [xLow(1:end-1) xLeft(1:end-1) xUp(1:end-1) xRight(1:end)];
    b(1:m+1) = [yLow(1:end-1) yLeft(1:end-1) yUp(1:end-1) yRight(1:end)];
    
    %rotation
    a(1:m+1) = (a(1:m+1)-l/2)*cos(theta) - (b(1:m+1)-R)*sin(theta) + l/2;
    b(1:m+1) = (a(1:m+1)-l/2)*sin(theta) + (b(1:m+1)-R)*cos(theta) + R;
    b(1:m+1) = b(1:m+1);
        
    %build the bubble interface at the first iteration, cirle
    if PARAM.alpha>PARAM.R
        
        %build the bubble interface at the first iteration
        PARAM.Vtot = 4/3*pi*R^3*PARAM.alpha^3;
        PARAM.h = PARAM.h + PARAM.thickness;
        [R1,R2,H] = findR1R2H(PARAM);
        [a(m+2:m+q+2),b(m+2:m+q+2)] = draw_test3(PARAM.start,0,H,R1,R2,q,theta);
        
    else
        
        [a(m+2:m+q+2),b(m+2:m+q+2)] = draw_circle_lean(0,0,PARAM.q,PARAM.D);
        a(m+2:m+q+2) = a(m+2:m+q+2)*PARAM.alpha + PARAM.start;
        b(m+2:m+q+2) = b(m+2:m+q+2)*PARAM.alpha;
        
    end
     
end




