%create two vectors containing the boundaries nodes (first point and last point of the element) once you have given the radius R, the
%lenght L of the pipe, the number of element in each region (n for the
%outlet, m for the wall, j for the inlet).

function [a,b] = droplet(R,q,R_bubble)

     if R_bubble>=R
         b = 0.9*R;
         a = R_bubble^3/b^2;
         D = (a-b)/(a+b);
     else
         D = 0;
     end
     
     %build the bubble interface 1 at the first iteration
     [a(1:q+1),b(1:q+1)] = draw_bubble(0,0,q,D,R_bubble);
     
%      figure
%      plot(a,b,'o-')
%      xlabel('x')
%      ylabel('r')
%      axis equal
%      title('domain')
%      hold on
end