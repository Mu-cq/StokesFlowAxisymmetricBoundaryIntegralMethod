%draw sphere pertubed by second legendre mode of a given volume

function [x,y] = spherePlus2LegendreCurvilinear(t,delta1,delta2,PARAM)

%legendre
P2 = legendre(2,cos(pi*t));
P2 = P2(1,:)';
P4 = legendre(4,cos(pi*t));
P4 = P4(1,:)';

%draw
xfun = @(f0) (f0+delta1*P2+delta2*P4).*cos(pi*t);
yfun = @(f0) (f0+delta1*P2+delta2*P4).*sin(pi*t);
V = VolumeCurvilinearAxisSpectral(xfun(1),xfun(1),PARAM);

%correct the volume with bisection
eps1 = 0.2; eps2 = 10;
V0 = PARAM.V0;
count = 0;
while abs(V0-V)>1e-15
    
   epsilon = (eps1+eps2)/2;
   x = xfun(epsilon);
   y = yfun(epsilon);
   V = VolumeCurvilinearAxisSpectral(x,y,PARAM);
   
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