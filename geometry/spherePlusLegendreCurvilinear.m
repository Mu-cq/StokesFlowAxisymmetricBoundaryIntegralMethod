%draw sphere pertubed by second legendre mode of a given volume

function [x,y] = spherePlusLegendreCurvilinear(t,delta,PARAM)

%legendre
P2 = legendre(2,cos(pi*t));
P2 = P2(1,:)';

%draw
x = (1+delta*P2).*cos(pi*t);
y = (1+delta*P2).*sin(pi*t);
V = VolumeCurvilinearAxisSpectral(x,y,PARAM);

%correct the volume with bisection
eps1 = 0.2; eps2 = 10;
V0 = PARAM.V0;
count = 0;
while abs(V0-V)>1e-15
    
   epsilon = (eps1+eps2)/2;
   x = epsilon*(1+delta*P2).*cos(pi*t);
   y = epsilon*(1+delta*P2).*sin(pi*t);
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