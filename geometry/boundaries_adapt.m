%create two vectors containing the boundaries nodes (first point and last point of the element) once you have given the radius R, the
%lenght L of the pipe, the number of element in each region (n for the
%outlet, m for the wall, j for the inlet).

function [a,b]=boundaries_adapt(PARAM)

    m = PARAM.m;    q = PARAM.q;
    R_bubble = PARAM.alpha*PARAM.R;
    L = PARAM.L;
    theta = PARAM.theta;
    R = PARAM.R;
    g = PARAM.tune_wall;

     if R_bubble >= PARAM.R
         b = 0.9*PARAM.R;
         a = R_bubble^3/b^2;
         D = (a-b)/(a+b);
     else
         D = 0;
     end         

     %build the panel at the wall STRAIGTH PART
     %[a(1:m+1),b(1:m+1)]=drawline2(PARAM.L,PARAM.R+PARAM.L/2*sin(PARAM.theta),0,PARAM.R-PARAM.L/2*sin(PARAM.theta),PARAM.m);
     s = 0:1/m:1;
     %adaptivity on one side
     %phi = atan(s)/atan(1); %this choose the mapping
     phi = (exp(-s).^g-1)/(exp(-1)^g-1); %this choose the mapping
     a(1:m+1) = L*(1 - phi);
     R1 = R + L/2*sin(theta);
     R2 = R - L/2*sin(theta);
     b(1:m+1) = (R2-R1)*(1 - phi) + R1;
     
     %build the bubble interface at the first iteration
     [a(m+2:m+q+2),b(m+2:m+q+2)] = draw_bubble(PARAM.start,0,PARAM.q,D,R_bubble);
     
%      figure
%      plot(a,b,'o-')
%      xlabel('x')
%      ylabel('r')
%      axis equal
%      title('domain')
%      hold on
end