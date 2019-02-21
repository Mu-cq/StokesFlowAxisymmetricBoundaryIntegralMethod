%compute green function in axisymmetric

function [SXX,SXY,SYX,SYY,Iaxis,Naxis] = MYaxisKernelSingleLayer (X,Y,X0,Y0)

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
 
  %end

 DX   = X-X0;
 DX2  = DX.*DX;
 DR2  = DX2+DY2;
 DR   = sqrt(DR2);
 
 % k parameter
 K2 = 4*Y.*Y0./(DX.^2+(Y+Y0).^2);
 K = sqrt(K2);
 K4 = K2.^2;

 [F, E] = ell_int_vect(K);
 %[F, E] = ellipke(K.*K);
 
 %prefactor
 alpha1 =2*K./sqrt(Y.*Y0);
 alpha3 = alpha1.^3;

 %compute I10
 I10 = alpha1.*F;
 I10(isnan(I10)) = 0;
 
 %compute I11
 num = (K2-2).*F - 2*E;
 I11 = alpha1.*num./K2;
 I11(isnan(I11)) = 0;
 
 %compute I30
 den = K.^2-1;
 I30 = alpha3.*E./den;
 I30(isnan(I30)) = 0;
 
 %compute I31
 num = 2*(1-K2).*F + (K2-2).*E;
 den = K2.*(K2-1);
 I31 = alpha3.*num./den;
 I31(isnan(I31)) = 0;
 
 %compute I32
 num = 4*(2-3*K2+K4).*F - (8-8*K2+K4).*E;
 den = K4.*(K2-1);
 I32 = alpha3.*num./den;
 I32(isnan(I32)) = 0;
   

 %  point x0 on the axis

 DR3  = DR2.*DR.*Iaxis;
 I10 = I10+pi2./DR.*Iaxis;
 temp = pi2./DR3;
 temp(isinf(temp)) = 0;
 I30 = I30+temp;
 temp = pi./DR3;
 temp(isinf(temp)) = 0;
 I32 = I32+temp;

    %end
 

%---
% build the Green's function
%---

 SXX = Y.*   (  I10 + DX2.*I30);
 SXY = Y.*DX.*(Y.*I30 - Y0 .*I31);
 SYX = Y.*DX.*(Y.*I31 - Y0 .*I30);
 SYY = Y.*   (  I11 + YY2.*I31-YYP.*(I30+I32));
 
%------
% also build the stress tensor
%------


%-----
% done
%-----

return
