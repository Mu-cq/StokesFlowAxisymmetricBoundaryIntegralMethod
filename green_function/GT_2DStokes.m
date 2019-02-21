%Stokeslet and Stresslet for 2D Stokes

function [SXX,SXY,SYY,QXXX,QXXY,QXYY,QYYY] ...
    = GT_2DStokes (X,Y,X0,Y0)

    %geometrical and often used variables
    r = sqrt((X-X0).^2+(Y-Y0).^2);
    dx = X-X0;
    dy = Y-Y0;
    
    %Stokeslet
    SXX = -log(r)+(dx.*dx)./r.^2;
    SXY = (dx.*dy)./r.^2;
    SYY = -log(r)+(dy.*dy)./r.^2;
    
    %Stresslet
    QXXX = -4*(dx.*dx.*dx)./r.^4;
    QXXY = -4*(dx.*dx.*dy)./r.^4;
    QXYY = -4*(dx.*dy.*dy)./r.^4;
    QYYY = -4*(dy.*dy.*dy)./r.^4;

end