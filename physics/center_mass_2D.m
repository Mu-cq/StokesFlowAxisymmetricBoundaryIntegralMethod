%compute center of mass of a 2D figure

function [xcm,ycm,A] = center_mass_2D(x,y)

        %normal vectors
        dx = x(2:end)-x(1:end-1);
        dy = y(2:end)-y(1:end-1);
        dl  = sqrt(dx.^2+dy.^2);
        n = [dy./dl -dx./dl];
        
        %area with divergence theorem (rectangle rule)
        int = [x/2 y/2];
        temp = sum((int(1:end-1,:)+int(2:end,:)).*n,2);
        A = sum(temp.*dl/2);
        
        %CM x coordinate (rectangle rule)
        int = [x.^2/2 zeros(numel(x),1)];
        temp = sum((int(1:end-1,:)+int(2:end,:)).*n,2);
        numx = sum(temp.*dl/2);
        
        xcm = numx/A;
        
        %CM y coordinate (rectangle rule)
        int = [zeros(numel(x),1) y.^2/2];
        temp = sum((int(1:end-1,:)+int(2:end,:)).*n,2);
        numy = sum(temp.*dl/2);
        
        ycm = numy/A;

return