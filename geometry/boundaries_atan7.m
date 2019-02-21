%create two vectors conteining the boundaries nodes (first point and last point of the element) once you have given the radius R, the
%lenght L of the pipe, the number of element in each region (n for the
%outlet, m for the wall, j for the inlet).

function [a,b]=boundaries_atan7(R1,R2,L,n,m,j,q)

     no = R1/R2;
     %no = 0.5;

     %build the panel at the outlet
     [a(1:ceil(n*no)+1),b(1:ceil(n*no)+1)]=drawline_atan(L,0,L,R1,ceil(n*no));
     
     %build the panel at the outlet
     [a(ceil(n*no)+1:n+1),b(ceil(n*no)+1:n+1)]=drawline_atan(L,R1,L,R2,n-ceil(n*no));

     %build the panel at the wall
     [a(n+1:m+n+1),b(n+1:m+n+1)]=drawline_atan(L,R2,0,R2,m);

     %build the panel at the inlet
     [a(n+m+1:n+m+ceil(j*(1-no))+1),b(n+m+1:n+m+ceil(j*(1-no))+1)]=drawline_atan(0,R2,0,R1,ceil(j*(1-no)));
     
     %build the panel at the inlet
     [a(n+m+ceil(j*(1-no))+1:n+m+j+1),b(n+m+ceil(j*(1-no))+1:n+m+j+1)]=drawline_atan(0,R1,0,0,j-ceil(j*(1-no)));

     %build the interface at the first iteration
     %[a(n+m+j+2:n+m+j+q+2),b(n+m+j+2:n+m+j+q+2)]=drawline2(0,R1,L,R1,q);
     [a(n+m+j+2:n+m+j+q+2),b(n+m+j+2:n+m+j+q+2)]=drawline_x3(0,R1,L,R1,q);
     
     figure
     plot(a,b,'o-')
     xlabel('x')
     ylabel('r')
     axis equal
     title('domain')
     hold on
end