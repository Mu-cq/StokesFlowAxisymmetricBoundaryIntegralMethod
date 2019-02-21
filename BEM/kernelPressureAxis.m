function [LX,LY,NXX,NXY,NYX,NYY] = kernelPressureAxis(X,Y,X0,Y0)

%LX and LY are the pressure in (X0,Y0) associated to a Stokeslet placed in (X,Y)

%NIXX, NIXY, NIYX, NIYY, are the pressure in (X0,Y0) associated with a
%Stresslet in (X,Y)

%eps=0.0000001;

pi2 = 2.0*pi;

 %on axis
 Iaxis = (Y0<eps);
   
 Y0(Y0<eps) = 0;

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
 
 %  point x0 off the axis
 
  YYR = sqrt(YYP);
  YY3 = YYR.^3;
  YY5 = YYR.^5;
  SY  = Y+Y0;
  SY2 = SY.^2;
 
  %end

 DX   = X-X0;
 %DX   = -X+X0;
 DX2  = DX.*DX;
 DR2  = DX2+DY2;
 DR   = sqrt(DR2);
 DXYY = DX2+YY2;
   
   %  point x0 off the axis

  FC1  = 4.0D0*YYP;
  FC2  = DX2+SY2;
  RK2  = FC1./FC2;
  RK   = sqrt(RK2);
  %RK3  = RK2.*RK;
  RK4  = RK2.^2;
  RK5  = RK4.*RK;
  RK2P = 1.0D0-RK2;

 [F, E] = ell_int_vect(RK);
 %[F, E] = ellipke(RK.*RK);

 %RJ10(isinf(RJ10)) = 0;
 RJ30 = 2.0D0*RK.*E./(YYR.*DR2);
 RJ30(isnan(RJ30)) = 0;
 RJ31 = RK.*(-F+DXYY.*E./DR2)./YY3;
 RJ31(isnan(RJ31)) = 0;

 

   %RL10 = F.*Naxis;
   RL30 = E./RK2P;
%  RK8  = RK6*RK2
   RL50 =  (2.0D0 * (2.0D0 - RK2).*RL30 - F) ./ (3.0D0*RK2P);
   RL52 =  (RL50 - RL30)./RK2;
   RL54 =  (F-RL50+2.0*RK2.*RL52)./RK4;
   
   FCTR = RK5./(8.0*YY5);
   FCTR(isinf(FCTR)) = 0;
   RJ50 = FCTR .* RL50;
   RJ50(isnan(RJ50)) = 0;
   RJ51 = FCTR .* (2.0D0 * RL52-RL50);
   RJ51(isnan(RJ51)) = 0;
   RJ52 = FCTR .* (4.0D0 * (RL54-RL52)+RL50);
   RJ52(isnan(RJ52)) = 0;
    %  point x0 on the axis

 DR3  = DR2.*DR.*Iaxis;
 DR5  = DR3.*DR2.*Iaxis;
 %RJ11 = 0.0;
 temp = pi2./DR3;
 temp(isinf(temp)) = 0;
 RJ30 = RJ30+temp;
 %RJ31 = 0.0;
 
 temp = pi2./DR5;
 temp(isinf(temp)) = 0;
 RJ50 = RJ50+temp;
 %RJ51 = 0.0;
 temp = pi./DR5;
 temp(isinf(temp)) = 0;
 RJ52 = RJ52+temp;

    %end
 

%Single layer
LX = -2*Y.*DX.*RJ30;
LY = -2*Y.*(Y.*RJ30-Y0.*RJ31);
 
%double layer
NXX = -4*Y.*(-RJ30 + 3*DX.^2.*RJ50);
NXY = -12*Y.*DX.*(Y.*RJ50 - Y0.*RJ51);
NYX = NXY;
NYY = -4*Y.*(-RJ30 + 3*Y0.^2.*RJ52 - 6*Y.*Y0.*RJ51 + 3*Y.^2.*RJ50);

%check I52 from Mathematica
%Int5 = ( -2*(4-6*RK.^2+RK.^6).*E+(8-16*RK.^2+7*RK.^4+RK.^6).*F  )./(3*RK.^4.*(-1+RK.^2).^2);
%INT5 = 4*RK.^5./(4.*Y.*Y0).^(2.5).*Int5;

return
