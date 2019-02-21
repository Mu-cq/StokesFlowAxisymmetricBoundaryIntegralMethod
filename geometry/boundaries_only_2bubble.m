%create two vectors containing the boundaries nodes (first point and last point of the element) once you have given the radius R, the
%lenght L of the pipe, the number of element in each region (n for the
%outlet, m for the wall, j for the inlet).

function [a,b]=boundaries_only_2bubble(R,q,p,R_bubble1,R_bubble2,dist)

     if R_bubble1>=R
         b = 0.9*R;
         a = R_bubble1^3/b^2;
         D1 = (a-b)/(a+b);
     else
         D1 = 0;
     end
     
     if R_bubble2>=R
         b = 0.9*R;
         a = R_bubble2^3/b^2;
         D2 = (a-b)/(a+b);
     else
         D2 = 0;
     end
     
     %build the bubble interface 1 at the first iteration
     [a(1:q+1),b(1:q+1)] = draw_bubble(dist/2,0,q,D1,R_bubble1);
     
     %build the bubble interface 2 at the first iteration
     [a(q+2:q+p+2),b(q+2:q+p+2)] = draw_bubble(-dist/2,0,p,D2,R_bubble2);
     
%      figure
%      plot(a,b,'o-')
%      xlabel('x')
%      ylabel('r')
%      axis equal
%      title('domain')
%      hold on
end