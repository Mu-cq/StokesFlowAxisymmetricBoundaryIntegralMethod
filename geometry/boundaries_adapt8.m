%create two vectors containing the boundaries nodes (first point and last point of the element) once you have given the radius R, the
%lenght L of the pipe, the number of element in each region (n for the
%outlet, m for the wall, j for the inlet).

function [a,b]=boundaries_adapt8(PARAM)

    m = PARAM.m;    q = PARAM.q;
    m1 = round(PARAM.LowUp*PARAM.m);     m2 = m-m1;    
    l = PARAM.L;
    theta = PARAM.theta;
    R = PARAM.R;
    g1low = PARAM.tune_wall;
    g1Up = 1;
    g2 = PARAM.tune_curve;
    thickness = PARAM.thickness;       

    %build lower and left part of the panel
    s = 0:1/m1:1;                                                % curvilinear coordinate
    s = (exp(-s).^g1low-1)/(exp(-1)^g1low-1)*(l+pi*thickness);   % this choose the mapping and arclenght
    x1 = 0;    y1 = 0;                                           % center of the small circle
    phi = (s-l)/thickness;                                       % angle for small circle
    aWallDown = (l-s).*(s<l) + (thickness.*cos(-phi+3*pi/2)+x1).*(s>=l);
    bWallDown = -thickness.*(s<l) + (thickness.*sin(-phi+3*pi/2)+y1).*(s>=l) + R;
     
    %build upper and right part of the panel
    s = (0:1/m2:1)*(l+pi*thickness);                       % curvilinear coordinate
    x1 = l;    y1 = 0;                                     % center of the small circle
    phi = (s-l)/thickness;                                 % angle for small circle
    aWallUp = s.*(s<l) + (thickness.*cos(-phi+pi/2)+x1).*(s>=l);
    bWallUp = thickness.*(s<l) + (thickness.*sin(-phi+pi/2)+y1).*(s>=l) + R;
     
    %paste the two wall
    a(1:m+1) = [aWallDown(1:end-1) aWallUp];
    b(1:m+1) = [bWallDown(1:end-1) bWallUp];
     
    %rotation
    a(1:m+1) = (a(1:m+1)-l/2)*cos(theta) - (b(1:m+1)-R)*sin(theta) + l/2;
    b(1:m+1) = (a(1:m+1)-l/2)*sin(theta) + (b(1:m+1)-R)*cos(theta) + R;
    b(1:m+1) = b(1:m+1);
        
    %build the bubble interface at the first iteration, cirle
    [a(m+2:m+q+2),b(m+2:m+q+2)] = draw_circle_lean(0,0,PARAM.q,PARAM.D);
    a(m+2:m+q+2) = a(m+2:m+q+2)*PARAM.alpha + PARAM.start;
    b(m+2:m+q+2) = b(m+2:m+q+2)*PARAM.alpha;
     
end