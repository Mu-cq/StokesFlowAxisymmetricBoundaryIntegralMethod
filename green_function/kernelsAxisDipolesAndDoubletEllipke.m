function [SXXD,SXYD,SYXD,SYYD,QXXXD,QXXYD,QXYXD,QXYYD,QYXXD,QYXYD,QYYXD,QYYYD,...
    SXXSD,SXYSD,SYXSD,SYYSD,QXXXSD,QXXYSD,QXYXSD,QXYYSD,QYXXSD,QYXYSD,QYYXSD,QYYYSD] ...
...
  = kernelsAxisDipolesAndDoubletEllipke(X,Y,X0,Y0)

error('BUG')

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
 %DY  = -Y+Y0;
 DY2 = DY.*DY;


  Y3   = Y.*Y2;
  Y66  = 6.0D0*Y;
  Y03  = Y0.*Y02;
 
 %  point x0 off the axis
 
  Y4  = Y2.*Y2;
  Y04 = Y02.*Y02;
  YYR = sqrt(YYP);
  YY3 = YYR.^3;
  YY5 = YYR.^5;
  SY  = Y+Y0;
  SY2 = SY.^2;
 
  %end

 DX   = X-X0;
 %DX   = -X+X0;
 DX2  = DX.*DX;
 DX4  = DX2.*DX2;
 DR2  = DX2+DY2;
 DR   = sqrt(DR2);
 DXYY = DX2+YY2;

 
 Y6DX  = Y66.*DX;
 Y6DX2 = Y66.*DX2;
 DX3   = DX.*DX2;
 Y6DX3 = Y66.*DX3;
   
   %  point x0 off the axis

  FC1  = 4.0D0*YYP;
  FC2  = DX2+SY2;
  RK2  = FC1./FC2;
  RK   = sqrt(RK2);
  %RK3  = RK2.*RK;
  RK4  = RK2.^2;
  RK5  = RK4.*RK;
  RK2P = 1.0D0-RK2;

 %[F, E] = ell_int_vect(RK);
 [F, E] = ellipke(RK.*RK);

 RJ10 = 2.0D0 * RK .* F./YYR;
 %RJ10(isinf(RJ10)) = 0;
 RJ10(isnan(RJ10)) = 0;
 RJ11 = RK.*(DXYY.*F-(DX2+SY2).*E)./YY3;
 RJ11(isnan(RJ11)) = 0;
 RJ30 = 2.0D0*RK.*E./(YYR.*DR2);
 RJ30(isnan(RJ30)) = 0;
 RJ31 = RK.*(-F+DXYY.*E./DR2)./YY3;
 RJ31(isnan(RJ31)) = 0;
 RJ32 = RK.*(-DXYY.*F+(DX4+2.0*DX2.*YY2+Y4+Y04).*E./DR2)./YY5;
 RJ32(isnan(RJ32)) = 0;

 

   %RL10 = F.*Naxis;
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
   FCTR(isinf(FCTR)) = 0;
   RJ50 = FCTR .* RL50;
   RJ50(isnan(RJ50)) = 0;
   RJ51 = FCTR .* (2.0D0 * RL52-RL50);
   RJ51(isnan(RJ51)) = 0;
   RJ52 = FCTR .* (4.0D0 * (RL54-RL52)+RL50);
   RJ52(isnan(RJ52)) = 0;
   RJ53 = FCTR .* (8.0D0 * RL56 - 12.0D0 *RL54 +6.0D0 * RL52 - RL50);
   RJ53(isnan(RJ53)) = 0;
%  RJ54 = FCTR * (16.0*RL58 - 32.0*RL56 + 24.0*RL54 ...
%                  -8.0*RL52 + RL50);

    %  point x0 on the axis

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
   %RJ53 = 0.0;
%  RJ54 = 3.0*pi/(4.0*DR5)

    %end
    
    RL70 =  (0.80*(2.0-RK2).*RL50 - 0.60*RL30) ./ RK2P;
    RL72 =  (RL70-RL50) ./ RK2;
    RL74 =  (RL30-RL70+2.0D0*RK2.*RL72) ./ RK4;
    RL76 = -(F-RL70+3.0D0*RK2.*RL72-3.0D0*RK4.*RL74) ./ RK6;
    RK7  = RK2.*RK5;
    YY7  = YYR.*YYR.*YY5;
    FCT7 = RK7./(32.0D0*YY7);
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
    %RJ71 = 0.0D0;
    temp = pi./DR7;
    temp(isinf(temp)) = 0;
    RJ72 = RJ72+temp;
    %RJ73 = 0.0D0;
 

%---
% build the Green's function
%---

 SXXD = Y.*(-RJ30+3*DX2.*RJ50);
 SXYD = 3*Y.*DX.*(Y.*RJ50-Y0.*RJ51);
 SYXD = 3*Y.*DX.*(Y.*RJ51 - Y0 .*RJ50);
 SYYD = Y.*(RJ31 -3*YY2.*RJ51-3*YYP.*(RJ50+RJ52));
 
 SXXSD = DX.*SXXD;
 SXYSD = DX.*SXYD - Y.*(Y.*RJ30-Y0.*RJ31);
 SYXSD = DX.*SYXD - Y.*(Y.*RJ31-Y0.*RJ30);
 SYYSD = DX.*SYYD;
 
%------
% also build the stress tensor
%------

  

  QXXXD = -Y66.*DX .*(-3*RJ50+5*DX2.*RJ70);
  QXXYD = Y66 .* (Y.*RJ50 - Y0.*RJ51-5*DX2.*(Y.*RJ70-Y0.*RJ71));
  QXYXD = QXXYD;
  QXYYD = Y6DX  .* (RJ50-5*(Y02.*RJ72+Y2.*RJ70 - 2.0*YYP.*RJ71));
  QYXXD = Y66.*(-Y.*RJ51+Y0.*RJ50+5*DX2.*(Y.*RJ71-Y0.*RJ70));
  QYXYD = Y6DX.*(5*(Y02+Y2).*RJ71 - 5*YYP.*(RJ72+RJ70)-RJ51);
  QYYXD = QYXYD;
  QYYYD = Y66 .* (-3*Y.*RJ51 + Y0.*(RJ50+2*RJ52) + 5*Y.*Y0.^2.*RJ73 + ...
                5*Y.*(Y.^2+2*Y0.^2).*RJ71 - 5*Y0.*(2*Y.^2+Y0.^2).*RJ72 - 5*Y.^2.*Y0.*RJ70);
              
  QXXXSD = DX.*QXXXD;
  QXXYSD =  DX.*QXXYD + 6*Y.*DX.*(Y.*RJ50 - Y0.*RJ51);
  QXYXSD = QXXYSD;
  QXYYSD = DX.*QXYYD - 6*Y.*(DX2.*RJ50 - Y0.^2.*RJ50 - Y0.^2.*RJ52 - 2*Y.*Y0.*RJ51);
  QYXXSD = DX.*QYXXD + 6*Y.*DX.*(Y.*RJ51 - Y0.*RJ50);
  QYXYSD = DX.*QYXYD;
  QYYXSD = QYXYSD;
  QYYYSD = DX.*QYYYD + 6*Y.*DX.*(Y.*RJ51 - Y0.*RJ50);

%-----
% done
%-----

return
