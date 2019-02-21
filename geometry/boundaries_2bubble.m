%create two vectors containing the boundaries nodes (first point and last point of the element) once you have given the radius R, the
%lenght L of the pipe, the number of element in each region (n for the
%outlet, m for the wall, j for the inlet).

function [a,b]=boundaries_2bubble(R,L,n,m,j,q,p,theta,R_bubble1,R_bubble2,dist)

     a1 = R_bubble1;
     a2 = R_bubble2;

     if R_bubble1>=R
         b1 = 0.9*R;
         a1 = R_bubble1^3/b1^2;
         D1 = (a1-b1)/(a1+b1);
     else
         D1 = 0;
     end
     
     if R_bubble2>=R
         b2 = 0.9*R;
         a2 = R_bubble2^3/b2^2;
         D2 = (a2-b2)/(a2+b2);
     else
         D2 = 0;
     end    

     %build the panel at the outlet
     [a(1:n+1),b(1:n+1)]=drawline2(L,0,L,R+L/2*sin(theta),n);

     %build the panel at the wall
     [a(n+1:m+n+1),b(n+1:m+n+1)]=drawline2(L,R+L/2*sin(theta),0,R-L/2*sin(theta),m);

     %build the panel at the inlet
     [a(n+m+1:n+m+j+1),b(n+m+1:n+m+j+1)]=drawline2(0,R-L/2*sin(theta),0,0,j);
     
     %build the bubble interface 1 at the first iteration
     [a(n+m+j+2:n+m+j+q+2),b(n+m+j+2:n+m+j+q+2)] = draw_bubble(0,0,q,D1,R_bubble1);
     
     %build the bubble interface 2 at the first iteration
     [a(n+m+j+q+3:n+m+j+q+p+3),b(n+m+j+q+3:n+m+j+q+p+3)] = draw_bubble(0,0,p,D2,R_bubble2);
     
     %move drop1
     a(n+m+j+2:n+m+j+q+2) = a(n+m+j+2:n+m+j+q+2) + a1 + dist + L/2;
     
     %move drop2
     a(n+m+j+q+3:n+m+j+q+p+3) = a(n+m+j+q+3:n+m+j+q+p+3) - a2 - dist + L/2;
     
%      figure
%      plot(a,b,'o-')
%      xlabel('x')
%      ylabel('r')
%      axis equal
%      title('domain')
%      hold on
end