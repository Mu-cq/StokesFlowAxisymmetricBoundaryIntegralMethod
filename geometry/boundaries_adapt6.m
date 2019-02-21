%create two vectors containing the boundaries nodes (first point and last point of the element) once you have given the radius R, the
%lenght L of the pipe, the number of element in each region (n for the
%outlet, m for the wall, j for the inlet).

function [a,b]=boundaries_adapt6(PARAM)

    m = PARAM.m;    q = PARAM.q;
    L = PARAM.L;
    theta = PARAM.theta;
    R = PARAM.R;
    g1 = PARAM.tune_wall;       

     %build the panel at the wall STRAIGTH PART
     %[a(1:m+1),b(1:m+1)]=drawline2(PARAM.L,PARAM.R+PARAM.L/2*sin(PARAM.theta),0,PARAM.R-PARAM.L/2*sin(PARAM.theta),PARAM.m);
     s = 0:1/m:1;
     %adaptivity on one side
     phi = (exp(-s).^g1-1)/(exp(-1)^g1-1); %this choose the mapping
     a(1:m+1) = L*(1 - phi);
     %R1 = R + L/2*sin(theta);
     %R2 = R - L/2*sin(theta);
     b(1:m+1) = R;
     
     %rotation
     a(1:m+1) = (a(1:m+1)-L/2)*cos(theta) - (b(1:m+1)-R)*sin(theta) + L/2;
     b(1:m+1) = (a(1:m+1)-L/2)*sin(theta) + (b(1:m+1)-R)*cos(theta) + R;
     
     if PARAM.alpha >= 1
     
        %build the bubble interface at the first iteration
        PARAM.Vtot = 4/3*pi*R^3*PARAM.alpha^3;
        [R1,R2,H] = findR1R2H(PARAM);
     
        [a(m+2:m+q+2),b(m+2:m+q+2)] = draw_test3(PARAM.start,0,H,R1,R2,q,theta);
     
     elseif PARAM.alpha < 1
        
        %build the bubble interface at the first iteration, cirle
        [a(m+2:m+q+2),b(m+2:m+q+2)] = draw_circle_lean(0,0,PARAM.q,0);
        a(m+2:m+q+2) = a(m+2:m+q+2)*PARAM.alpha + PARAM.start;
        b(m+2:m+q+2) = b(m+2:m+q+2)*PARAM.alpha;
         
     end
     
end