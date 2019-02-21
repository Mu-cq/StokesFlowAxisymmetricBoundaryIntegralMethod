%compute area contained by a closed with divergence theorem

function A = compute_area_2D(x,y)

        %normal vectors
        dx = x(2:end)-x(1:end-1);
        dy = y(2:end)-y(1:end-1);
        dl  = sqrt(dx.^2+dy.^2);
        n = [dy./dl -dx./dl];
        
        %area with divergence theorem (trapezi rule) 
        int = [x/2 y/2];
        temp = sum((int(1:end-1,:)+int(2:end,:)).*n,2);
        A = sum(temp.*dl/2);

return