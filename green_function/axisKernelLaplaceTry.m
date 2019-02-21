function [G,GX,GY,Iaxis,Naxis] = axisKernelLaplaceTry (X,Y,X0,Y0)

%-----------------------------------------
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

%-------------------------------------------------
%  Axisymmetric Green's function of Stokes flow
%  in free space
%
%  Let b be the strength of the point-force ring located at x;
%  then the induced velocity field is:
%
%  ux(x0) = Sxx(x,x0) * bx + Sxs(x,x0) * bs
%  us(x0) = Ssx(x,x0) * bx + Sss(x,x0) * bs
%
%  The kernel of the axisymmetric double-layer potential is:
%
%   Idlpx(x0) = ux * ( Qxxx * vnx + Qxxs * vns)
%             + us * ( Qxsx * vnx + Qxss * vns)
%
%   Idlps(x0) = ux * ( Qsxx * vnx + Qsxs * vns)
%             + us * ( Qssx * vnx + Qsss * vns)
%
%  arguments of Qxxx are (x,x0)
%
%  This is the flow due to a ring distribution of stresslets
%
%  Pij is used for the desingularization of the dlp
%
%------------
%
%  Iopt = 1 produces only the Green's function
%  Iopt = 2 produces the Green's function and the stress tensor
%-----------------------------------------------------


%----------
% constants
%----------

eps=0.0000001;

pi2 = 2.0*pi;

%----------
% initialize
%----------

%------------------------
% field point on x axis ?
%------------------------

 %Iaxis = 0;

 %on axis
   Iaxis = (Y0<eps);
   
   Y0(Y0<eps) = 0;
 
 %off axis
   Naxis = logical(1-Iaxis);

%-------
% launch
%-------

 Y2  = Y.*Y;
 Y02 = Y0.*Y0;
 YY2 = Y2+Y02;
 YYP = Y.*Y0;
 DY  = Y-Y0;
 DY2 = DY.*DY;
 
 YYR = sqrt(YYP);
 YY3 = YYR.^3;
 SY  = Y+Y0;
 SY2 = SY.^2;

 DX   = X-X0;
 DX2  = DX.*DX;
 DR2  = DX2+DY2;
 DR   = sqrt(DR2);
 DXYY = DX2+YY2;

 FC1  = 4.0D0*YYP;
 FC2  = DX2+SY2;
 RK2  = FC1./FC2;
 RK   = sqrt(RK2);

 %[F, E] = ell_int_vect(RK);
 [F, E] = ellipke(RK.*RK);
 
 RJ30 = 2.0D0*RK.*E./(YYR.*DR2);
 RJ30(isnan(RJ30)) = 0;
 RJ31 = RK.*(-F+DXYY.*E./DR2)./YY3;
 RJ31(isnan(RJ31)) = 0;

 DR3  = DR2.*DR.*Iaxis;
 temp = pi2./DR3;
 temp(isinf(temp)) = 0;
 RJ30 = RJ30+temp;
 
%---
% green's function
%---

 G = RK.*F./YYR/pi2;
 
%------
% gradient
%------

 GX = -DX.*RJ30/2/pi2;
 GY = -(Y.*RJ30-Y0.*RJ31)/2/pi2;

%-----
% done
%-----

return
