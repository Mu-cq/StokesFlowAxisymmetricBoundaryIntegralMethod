%create two vectors containing the boundaries nodes (first point and last point of the element) once you have given the radius R, the
%lenght L of the pipe, the number of element in each region (n for the
%outlet, m for the wall, j for the inlet).

function [a,b]=boundaries_adapt5(PARAM)

    %m1 = round(0.47*PARAM.m);     m2 = round(0.03*PARAM.m);    m = PARAM.m;    q = PARAM.q;
    %m3 = m1;    m4 = round(0.05*PARAM.m);
    m1 = round(0.47*PARAM.m);     m2 = round(0.5*PARAM.m)-m1;    m = PARAM.m;    q = PARAM.q;
    m3 = m1;    m4 = m2;
    L = PARAM.L;
    theta = PARAM.theta;
    R = PARAM.R;
    g1 = PARAM.tune_wall;
    g2 = PARAM.tune_curve;
    thickness = PARAM.thickness;       

     %build the panel at the wall STRAIGTH PART DOWN
     s = 0:1/m1:1;
     %adaptivity on one side
     phi = (exp(-s).^g1-1)/(exp(-1)^g1-1); %this choose the mapping
     a(1:m1+1) = L*(1 - phi);
     b(1:m1+1) = R;
     
     %build the panel at the wall CURVED PART LEFT
     s = 0:1/m2:1;
     %adaptivity on one side
     %phi = (exp(s).^g2-1)/(exp(1)^g2-1); %this choose the mapping
     phi = s;
     x1 = a(m1+1);
     y1 = b(m1+1) + thickness;
     a(m1+1:m1+m2+1) = thickness.*cos(-phi*pi+3*pi/2)+x1;
     b(m1+1:m1+m2+1) = thickness.*sin(-phi*pi+3*pi/2)+y1;
     
     %m3 = m1;   m4 = m2;
     %build the panel at the wall STRAIGTH PART UP
     s = 0:1/m3:1;
     %adaptivity on one side
     phi = (exp(s).^g1-1)/(exp(1)^g1-1); %this choose the mapping
     a(m1+m2+1:m1+m2+m3+1) = L*phi;
     b(m1+m2+1:m1+m2+m3+1) = R+thickness*2;
     
     %build the panel at the wall CURVED PART RIGHT
     s = 0:1/m4:1;
     phi = s;
     x1 = a(m1+m2+m3+1);
     y1 = b(m1+m2+m3+1) - thickness;
     a(m1+m2+m3+1:m1+m2+m3+m4+1) = thickness.*cos(-phi*pi+pi/2)+x1;
     b(m1+m2+m3+1:m1+m2+m3+m4+1) = thickness.*sin(-phi*pi+pi/2)+y1;
     
     %rotation
     a(1:m+1) = (a(1:m+1)-L/2)*cos(theta) - (b(1:m+1)-R)*sin(theta) + L/2;
     b(1:m+1) = (a(1:m+1)-L/2)*sin(theta) + (b(1:m+1)-R)*cos(theta) + R;
     b(1:m+1) = b(1:m+1)-thickness;
        
     %build the bubble interface at the first iteration, cirle
     [a(m+2:m+q+2),b(m+2:m+q+2)] = draw_circle_lean(0,0,PARAM.q,PARAM.D);
     a(m+2:m+q+2) = a(m+2:m+q+2)*PARAM.alpha + PARAM.start;
     b(m+2:m+q+2) = b(m+2:m+q+2)*PARAM.alpha;
     
end