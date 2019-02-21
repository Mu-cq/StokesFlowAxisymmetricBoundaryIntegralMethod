%draw sphere pertubed by second legendre mode of a given volume

function [x,y] = spherePlus3LegendreCartesian(theta,fff,whichFFF)

%coeff
f2 = fff(1);
f3 = fff(2);
f4 = fff(3);

%legendre
P2 = legendre(whichFFF(1),cos(theta));
P2 = P2(1,:)';
P3 = legendre(whichFFF(2),cos(theta));
P3 = P3(1,:)';
P4 = legendre(whichFFF(3),cos(theta));
P4 = P4(1,:)';

%draw
xfun = @(f0) (f0+f2*P2+f3*P3+f4*P4)'.*cos(theta);
yfun = @(f0) (f0+f2*P2+f3*P3+f4*P4)'.*sin(theta);
V = axis_int_gauss_vect(xfun(1),xfun(1));

%correct the volume with bisection
eps1 = 0.2; eps2 = 10;
V0 = 4/3*pi;
count = 0;
while abs(V0-V)>1e-15
    
   epsilon = (eps1+eps2)/2;
   x = xfun(epsilon);
   y = yfun(epsilon);
   V = axis_int_gauss_vect(x,y);
   
   if V>V0
       eps2 = epsilon;
   elseif V<V0
       eps1 = epsilon;
   end
   
   count = count+1;
   if count>1e5
        display(['The droplet is not initialized with the desired volume, err=' num2str(abs(V0-V))])
        break
   end
       
end