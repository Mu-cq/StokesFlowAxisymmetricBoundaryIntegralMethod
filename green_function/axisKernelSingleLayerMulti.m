function [SXX,SXY,SYX,SYY] = axisKernelSingleLayerMulti (X,Y,X0,Y0)

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

%-------
% launch
%-------

 Y2  = Y.*Y;
 Y02 = Y0.*Y0;
 YY2 = Y2+Y02;
 YYP = Y.*Y0;
 DY  = Y-Y0;
 DY2 = DY.*DY;
 
 %  point x0 off the axis
 
  Y4  = Y2.*Y2;
  Y04 = Y02.*Y02;
  YYR = sqrt(YYP);
  YY3 = YYR.*YYR.*YYR;
  YY5 = YYR.*YYR.*YYR.*YYR.*YYR;
  SY  = Y+Y0;
  SY2 = SY.^2;
 
  %end

 DX   = X-X0;
 DX2  = DX.*DX;
 DX4  = DX2.*DX2;
 DR2  = DX2+DY2;
 DXYY = DX2+YY2;
   
   %  point x0 off the axis

  FC1  = 4*YYP;
  FC2  = DX2+SY2;
  RK2  = FC1./FC2;
  RK   = sqrt(RK2);

 %[F, E] = ell_int_vect(RK);
 [F, E] = ellipke(RK.*RK);

 RJ10 = 2 * RK .* F./YYR;
 RJ11 = RK.*(DXYY.*F-(DX2+SY2).*E)./YY3;
 RJ30 = 2*RK.*E./(YYR.*DR2);
 RJ31 = RK.*(-F+DXYY.*E./DR2)./YY3;
 RJ32 = RK.*(-DXYY.*F+(DX4+2.0*DX2.*YY2+Y4+Y04).*E./DR2)./YY5;

%---
% build the Green's function
%---

 SXX = Y.*   (  RJ10 + DX2.*RJ30);
 SXY = Y.*DX.*(Y.*RJ30 - Y0 .*RJ31);
 SYX = Y.*DX.*(Y.*RJ31 - Y0 .*RJ30);
 SYY = Y.*   (  RJ11 + YY2.*RJ31-YYP.*(RJ30+RJ32));
 
%------
% also build the stress tensor
%------


%-----
% done
%-----
