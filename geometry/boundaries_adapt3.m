%create two vectors containing the boundaries nodes (first point and last point of the element) once you have given the radius R, the
%lenght L of the pipe, the number of element in each region (n for the
%outlet, m for the wall, j for the inlet).

function [a,b]=boundaries_adapt3(PARAM)

    m1= round(0.8*PARAM.m);     m2 = PARAM.m-m1;    m = PARAM.m;    q = PARAM.q;
    R_bubble = PARAM.alpha*PARAM.R;
    L = PARAM.L;
    theta = PARAM.theta;
    R = PARAM.R;
    g1 = PARAM.tune_wall;
    g2 = PARAM.tune_curve;

     if R_bubble >= PARAM.R
         b = 0.9*PARAM.R;
         a = R_bubble^3/b^2;
         D = (a-b)/(a+b);
     else
         D = 0;
     end         

     %build the panel at the wall STRAIGTH PART
     %[a(1:m+1),b(1:m+1)]=drawline2(PARAM.L,PARAM.R+PARAM.L/2*sin(PARAM.theta),0,PARAM.R-PARAM.L/2*sin(PARAM.theta),PARAM.m);
     s = 0:1/m1:1;
     %adaptivity on one side
     %phi = atan(s)/atan(1); %this choose the mapping
     phi = (exp(-s).^g1-1)/(exp(-1)^g1-1); %this choose the mapping
     a(1:m1+1) = L*(1 - phi);
     R1 = R + L/2*sin(theta);
     R2 = R - L/2*sin(theta);
     b(1:m1+1) = (R2-R1)*(1 - phi) + R1;
     
     %build the panel at the wall CURVED PART
     %[a(1:m+1),b(1:m+1)]=drawline2(PARAM.L,PARAM.R+PARAM.L/2*sin(PARAM.theta),0,PARAM.R-PARAM.L/2*sin(PARAM.theta),PARAM.m);
     s = 0:1/m2:1;
     %adaptivity on one side
     %phi = atan(s)/atan(1); %this choose the mapping
     phi = (exp(s).^g2-1)/(exp(1)^g2-1); %this choose the mapping
     x1 = a(m1+1);
     y1 = b(m1+1) + R/10;
     a(m1+1:m1+m2+1) = R/10.*cos(-phi*pi+3*pi/2)+x1;
     b(m1+1:m1+m2+1) = R/10.*sin(-phi*pi+3*pi/2)+y1;
     
     %build the bubble interface at the first iteration
     h = PARAM.h;
     Rb = PARAM.alpha*R;
     r = R-h;
     Vr = 4/3*pi*(Rb^3-r^3);
     H = Vr/pi/r^2;
     [a(m+2:m+q+2),b(m+2:m+q+2)] = draw_test2(PARAM.start,0,H,r,q);
     
%      figure
%      plot(a,b,'o-')
%      xlabel('x')
%      ylabel('r')
%      axis equal
%      title('domain')
%      hold on
end