function [SXX,SXY,SYX,SYY ...
         ,QXXX,QXXY,QXYX,QXYY ...
         ,QYXX,QYXY,QYYX,QYYY] = sgf_ax_w_vect (X,Y,X0,Y0,wall)

%-----------------------------------------
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

%---------------------------------------------------
%  Axisymmetric Green's function of Stokes flow
%  and kernel of the double-layer potential
%  for flow bounded by a plane wall
%  located at x = wall
%
%  (x,   y) are the x and s coordinates of the singularity
%  (x0, y0) are the x and s coordinates of the field point
%
%  Sij     Green's function
%  Qijk    Double-layer kernel
%
%  Iaxis = 0  field point is off-axis s = 0
%        = 1  field point is on-axis
%
%  wall:   x location of the wall
%
%  h:  distance of field point from wall
%
%  Use with the boundary integral equation:
%  ---------------------------------------
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
%  The arguments of Qxxx are (x,x0)
%  This is the flow due to a ring distribution of stresslets
%
%  Pij is used for the desingularization of the dlp
%
%------------
%  Iopt = 1 produces only the Green's function
%  Iopt = 2 produces the Green's function and the stress tensor
%-------------------------------------------------

%----------
% initialize
%----------

QXXX=0; QXXY=0; QXYX=0;QXYY=0;
QYXX=0; QYXY=0; QYYX=0;QYYY=0;
PXX=0;  PXY=0;  PYX=0;PYY=0;

%----------
% constants
%----------

 eps=0.00000001;

 pi2 = 2.0D0*pi;

%----------------------
% field point on axis ?
%----------------------

 Iaxis = (Y0<eps);
 Y0(Y0<eps) = 0;
 
 %off axis
 Naxis = logical(1-Iaxis);

%---------------
% prepare to run
%---------------

  H    = X0-wall;
  HH   = 2.0D0*H;
  H2   = H.*H;
  H2H2 = H2+H2;

  Y2  = Y.*Y;
  Y02 = Y0.*Y0;
  YY2 = Y2+Y02;
  YYP = Y.*Y0;
  DY  = Y-Y0;
  DY2 = DY.*DY;

  Y3  = Y.*Y2;
  Y66 = 6.0D0*Y;
  Y03 = Y0.*Y02;
 
  Y4  = Y2.*Y2;
  Y04 = Y02.*Y02;
  YYR = sqrt(YYP);
  YY3 = YYR.*YYR.*YYR;
  YY5 = YYR.*YYR.*YY3;
  SY  = Y+Y0;
  SY2 = SY.*SY;
 
  DX = X-X0;

%-----------------------------
% will pass twice:
%
%  first for the primary ring
%  second for the image ring
%-----------------------------

%---
 for Ipass=1:2
%---

      DX2  = DX.*DX;
      DX4  = DX2.*DX2;
      DR2  = DX2+DY2;
      DR   = sqrt(DR2);
      DXYY = DX2+YY2;

      Y6DX  = Y66.*DX;
      Y6DX2 = Y66.*DX2;
      DX3   = DX.*DX2;
      Y6DX3 = Y66.*DX3;

 
%----------------------

        FC1  = 4.0D0*YYP;
        FC2  = DX2+SY2;
        RK2  = FC1./FC2;
        RK   = sqrt(RK2);
        RK4  = RK2.*RK2;
        RK5  = RK4.*RK;
        RK2P = 1.0D0-RK2;

        [F, E]=ell_int_vect(RK);

        RJ10 = 2.0D0 * RK .* F./YYR;
        RJ10(isnan(RJ10)) = 0;
        RJ11 = RK.*(DXYY.*F-(DX2+SY2).*E)./YY3;
        RJ11(isnan(RJ11)) = 0;
        RJ30 = 2.0D0*RK.*E./(YYR.*DR2);
        RJ30(isnan(RJ30)) = 0;
        RJ31 = RK.*(-F+DXYY.*E./DR2)./YY3;
        RJ31(isnan(RJ31)) = 0;
        RJ32 = RK.*(-DXYY.*F+(DX4+2.0*DX2.*YY2+Y4+Y04).*E./DR2)./YY5;
        RJ32(isnan(RJ32)) = 0;

        RL30 = E./RK2P;
        RK6  = RK4.*RK2;
    %   RK8  = RK6*RK2
        RL50 =  (2.0D0 * (2.0D0 - RK2).*RL30 - F) ./ (3.0D0*RK2P);
        RL52 =  (RL50 - RL30)./RK2;
        RL54 =  (F-RL50+2.0*RK2.*RL52)./RK4;
        RL56 = -(E-RL50+3.0*RK2.*RL52-3.0*RK4.*RL54)./ RK6;
    %   PREP = ( 2.0*(2.0-RK2)*E - RK2P*RK2P*F ) / 3.0D0
    %   RL58 = - (PREP - RL50 + 4.0D0 * RK2*RL52 - 6.0D0 * RK4*RL54 ...
    %                          + 4.0D0 * RK6*RL56 ) / RK8;
        FCTR = RK5./(8.0*YY5);
        FCTR(isinf(FCTR)) = 0;
        RJ50 = FCTR .* RL50;
        RJ50(isnan(RJ50)) = 0;
        RJ51 = FCTR .* (2.0D0 * RL52-RL50);
        RJ51(isnan(RJ51)) = 0;
        RJ52 = FCTR .* (4.0D0 * (RL54-RL52)+RL50);
        RJ52(isnan(RJ52)) = 0;
        RJ53 = FCTR .* (8.0D0 * RL56 - 12.0D0 *RL54 +6.0D0 * RL52 - RL50);
        RJ53(isnan(RJ53)) = 0;


%---
%---

         DR3  = DR2.*DR.*Iaxis;
         DR5  = DR3.*DR2.*Iaxis;
         RJ10 = RJ10+pi2./DR.*Iaxis;
         %RJ11 = 0.0;
         temp = pi2./DR3;
         temp(isinf(temp)) = 0;
         RJ30 = RJ30+temp;
         %RJ31 = 0.0;
         temp = pi./DR3;
         temp(isinf(temp)) = 0;
         RJ32 = RJ32+temp;

       temp = pi2./DR5;
       temp(isinf(temp)) = 0;
       RJ50 = RJ50+temp;
       %RJ51 = 0.0;
       temp = pi./DR5;
       temp(isinf(temp)) = 0;
       RJ52 = RJ52+temp;

%---
%---

%-----------------
% Green's function
%-----------------
 
      SXX = Y.*    (  RJ10+DX2.*RJ30);
      SXY = Y.* DX.*(Y.*RJ30-Y0 .*RJ31);
      SYX = Y.* DX.*(Y.*RJ31-Y0 .*RJ30);
      SYY = Y.*    (  RJ11+YY2.*RJ31-YYP.*(RJ30+RJ32));


        QXXX = - Y6DX3 .* RJ50;
        QXXY = - Y6DX2 .* (Y.*RJ50 - Y0.*RJ51);
        QXYX =   QXXY;
        QXYY = - Y6DX  .* (Y02.*RJ52+Y2.*RJ50 - 2.0*YYP.*RJ51);
        QYXX = - Y6DX2 .* (Y.*RJ51-Y0.*RJ50);
        QYXY = - Y6DX  .* ((Y02+Y2).*RJ51 - YYP.*(RJ52+RJ50));
        QYYX = QYXY;
        QYYY = - Y66 .* ( Y3     .*  RJ51 ... 
                       - Y03    .*  RJ52  ...
                       - Y.*YYP  .* (RJ50+2.0*RJ52) ...
                       + Y0.*YYP .* (RJ53+2.0*RJ51) );

%------
% build
%------

      if(Ipass==1)
        SXXSAVE = SXX;
        SXYSAVE = SXY;
        SYXSAVE = SYX;
        SYYSAVE = SYY;
        QXXXSAVE = QXXX;
        QXXYSAVE = QXXY;
        QXYXSAVE = QXYX;
        QYXXSAVE = QYXX;
        QXYYSAVE = QXYY;
        QYXYSAVE = QYXY;
        QYYXSAVE = QYYX;
        QYYYSAVE = QYYY;
        DX = X+X0-2.0D0*wall;
      else
        SXX = SXXSAVE - SXX;
        SXY = SXYSAVE - SXY;
        SYX = SYXSAVE - SYX;
        SYY = SYYSAVE - SYY;
        QXXX = QXXXSAVE - QXXX;
        QXXY = QXXYSAVE - QXXY;
        QXYX = QXYXSAVE - QXYX;
        QYXX = QYXXSAVE - QYXX;
        QXYY = QXYYSAVE - QXYY;
        QYXY = QYXYSAVE - QYXY;
        QYYX = QYYXSAVE - QYYX;
        QYYY = QYYYSAVE - QYYY;
      end

%---
end    % of Ipass
%---
 
%------------------------------------
% source dipole and stokeslet doublet
%------------------------------------

      DDXX = Y .*        (-RJ30+3.0*DX2.*RJ50);
      DDXY = Y.*3.0.*DX.*( Y.*RJ50-Y0.*RJ51);
      DDYX = Y .* 3.0.*DX.*(-Y.*RJ51+Y0.*RJ50);
      DDYY = Y .*     (RJ31-3.0*YY2.*RJ51+3.0*YYP.*(RJ50+RJ52));
 
      EEXX = DX.*DDXX;
      EEXY = DX.*DDXY + Y.*(Y0.*RJ31-Y.*RJ30);
      EEYX = DX.*DDYX + Y.*(Y0.*RJ30-Y.*RJ31);
      EEYY = DX.*DDYY;
 
      SXX = SXX + H2H2.*DDXX - HH.*EEXX;
      SXY = SXY + H2H2.*DDXY - HH.*EEXY;
      SYX = SYX + H2H2.*DDYX - HH.*EEYX;
      SYY = SYY + H2H2.*DDYY - HH.*EEYY;

%--- 
%---

        RL70 =  (0.80*(2.0-RK2).*RL50 - 0.60*RL30) ./ RK2P;
        RL72 =  (RL70-RL50) ./ RK2;
        RL74 =  (RL30-RL70+2.0D0*RK2.*RL72) ./ RK4;
        RL76 = -(F-RL70+3.0D0*RK2.*RL72-3.0D0*RK4.*RL74) ./ RK6;
        RK7  = RK2.*RK5;
        YY7  = YYR.*YYR.*YY5;
        FCT7 = RK7./(32.0D0*YY7);
        FCT7(isinf(FCT7)) = 0;
        RJ70 = FCT7.* RL70;
        RJ70(isnan(RJ70)) = 0;
        RJ71 = FCT7.*(2.0D0*RL72-RL70);
        RJ71(isnan(RJ71)) = 0;
        RJ72 = FCT7.*(4.0D0*(RL74-RL72)+RL70);
        RJ72(isnan(RJ72)) = 0;
        RJ73 = FCT7.*(8.0D0*RL76-12.0D0*RL74+6.0D0*RL72-RL70);
        RJ73(isnan(RJ73)) = 0;
      
        DR7  = DR5.*DR2.*Iaxis;
        temp = pi2./DR7;
        temp(isinf(temp)) = 0;
        RJ70 = RJ70+temp;
        temp = pi./DR7;
        temp(isinf(temp)) = 0;
        RJ72 = RJ72+temp;
      
%---
 
      DQXXX = - Y6DX .* (-3.0*RJ50 + 5.0*DX2.*RJ70);
      DQXXY = - Y66  .* (-Y.*RJ50+Y0.*RJ51 - 5.0*DX2.*(-Y.*RJ70+Y0.*RJ71));
      DQXYX =   DQXXY;
      DQXYY = - Y6DX .* (-RJ50 + 5.0*(Y2.*RJ70-2.0*YYP.*RJ71+Y02.*RJ72));
      DQYXX =   Y66  .* (-Y.*RJ51+Y0.*RJ50 - 5.0*DX2.*(-Y.*RJ71+Y0.*RJ70));
      DQYXY =   Y6DX .* (-RJ51 + 5.0*YY2.*RJ71 ...
                              - 5.0*YYP.*(RJ70+RJ72));
      DQYYX = DQYXY;
      DQYYY =    Y66.*(-3.0D0*Y.*RJ51+Y0.*(RJ50+2.0D0*RJ52) ...
                      +5.0D0*Y.*Y02.*RJ73 ...
                      +5.0D0*Y.*(Y2+2.0*Y02).*RJ71 ...
                      -5.0D0*Y0.*(2.0*Y2+Y02).*RJ72 ...
                      -5.0D0*Y2.*Y0.*RJ70);
 
      FQXXX =   0.0D0;
      FQXXY =   Y6DX .* (Y.*RJ50 - Y0.*RJ51);
      FQXYX =   FQXXY;
      FQXYY = - Y66  .*(DX2.*RJ50-Y2.*RJ50-Y02.*RJ52+2.0*YYP.*RJ51);
      FQYXX =   Y6DX .*(Y.*RJ51 - Y0.*RJ50);
      FQYXY =   0.0D0;
      FQYYX =   0.0D0;
      FQYYY =   Y6DX .*(Y.*RJ51 - Y0.*RJ50);
 
      EQXXX = DX.*DQXXX + FQXXX;
      EQXXY = DX.*DQXXY + FQXXY;
      EQXYX = DX.*DQXYX + FQXYX;
      EQYXX = DX.*DQYXX + FQYXX;
      EQXYY = DX.*DQXYY + FQXYY;
      EQYXY = DX.*DQYXY + FQYXY;
      EQYYX = DX.*DQYYX + FQYYX;
      EQYYY = DX.*DQYYY + FQYYY;
 
      QXXX = QXXX + H2H2.*DQXXX - HH.*EQXXX;
      QXXY = QXXY + H2H2.*DQXXY - HH.*EQXXY;
      QXYX = QXYX + H2H2.*DQXYX - HH.*EQXYX;
      QYXX = QYXX + H2H2.*DQYXX - HH.*EQYXX;
      QXYY = QXYY + H2H2.*DQXYY - HH.*EQXYY;
      QYXY = QYXY + H2H2.*DQYXY - HH.*EQYXY;
      QYYX = QYYX + H2H2.*DQYYX - HH.*EQYYX;
      QYYY = QYYY + H2H2.*DQYYY - HH.*EQYYY;

%---
 
%-----
% done
%-----

  return
