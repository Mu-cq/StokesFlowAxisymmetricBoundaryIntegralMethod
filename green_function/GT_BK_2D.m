% compute Stokeslet and Stresslet for Brinkman equation

function [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY] ...
    = GT_BK_2D (X,Y,X0,Y0,k)

    %geometrical and often used variables
    r = sqrt((X-X0).^2+(Y-Y0).^2);
    arg = k*r;
    dx = X-X0;
    dy = Y-Y0;

    %Bessel function of second kind
    K0 = besselk(0,arg);
    K1 = besselk(1,arg);

    %often used terms
    A1 = 2*(K0+K1./arg-1./arg.^2);
    A2 = 2*(-K0-2*K1./arg+2./arg.^2);
    
    %Stokeslet
    SXX = 1/4/pi*(A1+A2.*dx.*dx./r.^2);
    SXY = 1/4/pi*A2.*dx.*dy./r.^2;
    SYX = SXY;
    SYY = 1/4/pi*(A1+A2.*dy.*dy./r.^2);
    
    %Stresslet
    QXXX = dx.*(A2-1)/2/pi./r.^2+dx/pi./r.^2.*(A2-K1.*arg)-dx.^3/pi./r.^4.*(2*A2-K1.*arg);
    QXXY = dx/2/pi./r.^2.*(A2-K1.*arg)-dx.^2.*dy/pi./r.^4.*(2*A2-K1.*arg);
    QXYX = QXXY;
    QXYY = dx.*(A2-1)/2/pi./r.^2-dx.*dy.^2/pi./r.^4.*(2*A2-K1.*arg);
    QYXX = dy.*(A2-1)/2/pi./r.^2-dx.^2.*dy/pi./r.^4.*(2*A2-K1.*arg);
    QYXY = dx/2/pi./r.^2.*(A2-K1.*arg)-dx.*dy.^2/pi./r.^4.*(2*A2-K1.*arg);
    QYYX = QYXY;
    QYYY = dy.*(A2-1)/2/pi./r.^2+dy/pi./r.^2.*(A2-K1.*arg)-dy.^3/pi./r.^4.*(2*A2-K1.*arg);

end