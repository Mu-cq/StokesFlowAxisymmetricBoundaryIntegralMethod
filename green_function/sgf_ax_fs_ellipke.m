function [SXX,SXY,SYX,SYY ...
...
   ,QXXX,QXXY,QXYX,QXYY ...
   ,QYXX,QYXY,QYYX,QYYY ...
...
   ,PXX,PXY,PYX,PYY ...
   ,Iaxis] ...
...
  = sgf_ax_fs_ellipke (Iopt,X,Y,X0,Y0)

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

  QXXX = zeros(numel(X),1);
  QXXY = zeros(numel(X),1);
  QXYX = zeros(numel(X),1);
  QXYY = zeros(numel(X),1);
  QYXX = zeros(numel(X),1);
  QYXY = zeros(numel(X),1);
  QYYX = zeros(numel(X),1);
  QYYY = zeros(numel(X),1);

  PXX = zeros(numel(X),1);
  PXY = zeros(numel(X),1);
  PYX = zeros(numel(X),1);
  PYY = zeros(numel(X),1);

%------------------------
% field point on x axis ?
%------------------------

 Iaxis = 0;

 if(Y0<eps) 
   Iaxis = 1;
 end

%-------
% launch
%-------

 Y2  = Y.*Y;
 Y02 = Y0.*Y0;
 YY2 = Y2+Y02;
 YYP = Y.*Y0;
 DY  = Y-Y0;
 DY2 = DY.*DY;

 if(Iopt > 1)
  Y3   = Y.*Y2;
  Y66  = 6.0D0*Y;
  Y03  = Y0.*Y02;
 end
 
 if(Iaxis==0) % off the axis
  Y4  = Y2.*Y2;
  Y04 = Y02.*Y02;
  YYR = sqrt(YYP);
  YY3 = YYR.^3;
  YY5 = YYR.^5;
  SY  = Y+Y0;
  SY2 = SY.^2;
 end

 DX   = X-X0;
 DX2  = DX.*DX;
 DX4  = DX2.*DX2;
 DR2  = DX2+DY2;
 DR   = sqrt(DR2);
 DXYY = DX2+YY2;

 if(Iopt > 1)
   Y6DX  = Y66.*DX;
   Y6DX2 = Y66.*DX2;
   DX3   = DX.*DX2;
   Y6DX3 = Y66.*DX3;
 end

%------------------------
  if(Iaxis==0)  %  point x0 off the axis
%------------------------

  FC1  = 4.0D0*YYP;
  FC2  = DX2+SY2;
  RK2  = FC1./FC2;
  RK   = sqrt(RK2);
  RK3  = RK2.*RK;
  RK4  = RK2.^2;
  RK5  = RK4.*RK;
  RK2P = 1.0D0-RK2;

 [F, E] = ellipke(RK.*RK);

 RJ10 = 2.0D0 * RK .* F./YYR;
 RJ11 = RK.*(DXYY.*F-(DX2+SY2).*E)./YY3;
 RJ30 = 2.0D0*RK.*E./(YYR.*DR2);
 RJ31 = RK.*(-F+DXYY.*E./DR2)./YY3;
 RJ32 = RK.*(-DXYY.*F+(DX4+2.0*DX2.*YY2+Y4+Y04).*E./DR2)./YY5;

 if(Iopt ~=1)

   RL10 = F;
   RL30 = E./RK2P;
   RK6  = RK4.*RK2;
%  RK8  = RK6*RK2
   RL50 =  (2.0D0 * (2.0D0 - RK2).*RL30 - F) ./ (3.0D0*RK2P);
   RL52 =  (RL50 - RL30)./RK2;
   RL54 =  (F-RL50+2.0*RK2.*RL52)./RK4;
   RL56 = -(E-RL50+3.0*RK2.*RL52-3.0*RK4.*RL54)./ RK6;
%  PREP = ( 2.0*(2.0-RK2)*E - RK2P*RK2P*F ) / 3.0D0
%  RL58 = - (PREP - RL50 + 4.0D0 * RK2*RL52 - 6.0D0 * RK4*RL54 ...
%                          + 4.0D0 * RK6*RL56 ) / RK8;
   FCTR = RK5./(8.0*YY5);
   RJ50 = FCTR .* RL50;
   RJ51 = FCTR .* (2.0D0 * RL52-RL50);
   RJ52 = FCTR .* (4.0D0 * (RL54-RL52)+RL50);
   RJ53 = FCTR .* (8.0D0 * RL56 - 12.0D0 *RL54 +6.0D0 * RL52 - RL50);
%  RJ54 = FCTR * (16.0*RL58 - 32.0*RL56 + 24.0*RL54 ...
%                  -8.0*RL52 + RL50);

end

%-----------
   else       %  point x0 on the axis
%-----------

 DR3  = DR2.*DR;
 DR5  = DR3.*DR2;
 RJ10 = pi2/DR;
 RJ11 = 0.0;
 RJ30 = pi2/DR3;
 RJ31 = 0.0;
 RJ32 = pi/DR3;

 if(Iopt > 1)
   RJ50 = pi2/DR5;
   RJ51 = 0.0;
   RJ52 = pi/DR5;
   RJ53 = 0.0;
%  RJ54 = 3.0*pi/(4.0*DR5)
 end

%-------------
  end
%-------------

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

if(Iopt>1)

  QXXX = - Y6DX3 .* RJ50;
  QXXY = - Y6DX2 .* (Y.*RJ50 - Y0.*RJ51);
  QXYX =   QXXY;
  QXYY = - Y6DX  .* (Y02.*RJ52+Y2.*RJ50 - 2.0*YYP.*RJ51);
  QYXX = - Y6DX2 .* (Y  .*RJ51-Y0.*RJ50);
  QYXY = - Y6DX  .* ((Y02+Y2).*RJ51 - YYP.*(RJ52+RJ50));
  QYYX = QYXY;
  QYYY = - Y66 .* ( Y3     .*  RJ51 ...
                  - Y03    .*  RJ52 ...
                  - Y.*YYP  .* (RJ50+2.0*RJ52)...
                  + Y0.*YYP .* (RJ53+2.0*RJ51) );
  PXX = QYXX;
  PXY = QYXY;
  PYX = - Y6DX .* (Y2.*RJ52 - 2.0D0 * YYP.*RJ51 + Y02.*RJ50);
  PYY = - Y66  .* (  Y0.*YYP .* (2.0D0 * RJ52+RJ50) ...
                    - Y03    .* RJ51 ...
                    - Y .*YYP .* (2.0*RJ51+RJ53) ...
                    + Y3     .* RJ52  );
end

%-----
% done
%-----

return
