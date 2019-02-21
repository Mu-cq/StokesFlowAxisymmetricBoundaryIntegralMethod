%create two vectors conteining the boundaries nodes (first point and last point of the element) once you have given the radius R, the
%lenght L of the pipe, the number of element in each region (n for the
%outlet, m for the wall, j for the inlet).

function [a,b]=boundaries_atan(R,L,n,m,j)

     %build the panel at the outlet
     [a(1:n+1),b(1:n+1)]=drawline_atan(L,0,L,R,n);

     %build the panel at the wall
     [a(n+1:m+n+1),b(n+1:m+n+1)]=drawline_atan(L,R,0,R,m);

     %build the panel at the inlet
     [a(n+m+1:n+m+j+1),b(n+m+1:n+m+j+1)]=drawline_atan(0,R,0,0,j);
     
     figure
     plot(a,b)
     axis equal
     xlabel('x')
     ylabel('r')
     title('domain')
     hold on
end