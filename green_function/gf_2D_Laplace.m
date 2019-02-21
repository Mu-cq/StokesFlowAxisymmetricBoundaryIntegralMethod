%compute Green's and Neumann function for Laplace Laplace equation

function [G,Gx,Gy,NXX,NXY,NYX,NYY,Gxx] = gf_2D_Laplace(x,y,x0,y0,wall1,wall2,visc,p_out,mass_flow)

    %I WANT A SOURCE WHEN mass_flow IS POSITIVE
    mass_flow = -mass_flow;

%constants

      pi2 = 2*pi;
      pi4 = 4*pi;

%prepare

      h  = wall2-wall1;
      h2 = 2*h;

      wn = pi2/h2;        %wave number

      %sign = 1;                      %Green
%       If(Ign.eq.2) sign = -1.0D0        %Neumann

%primary array 

      dx = x-x0;
      dy = y-y0;

      A = wn*dx;
      B = wn*dy;
      C = cosh(A)-cos(B);

%image system
      
      y0i = 2.0*wall1-y0;

      dxi = dx;
      dyi = y-y0i;

      Ai = wn*dxi;
      Bi = wn*dyi;
      Ci = cosh(Ai)-cos(Bi);

      G = ( -log(C) - log(Ci) )/pi4;
      
      G = G*mass_flow;

%compute the gradient

%       cf = -1/(2 * h2);
% 
%       %negative mass injection (I modified the signs)
%       Gx = cf * (sinh(A)./C + sinh(Ai)./Ci);
%       Gy = cf * (sin (B)./C +  sin(Bi)./Ci);

        %MATHEMATICA
        Gx = (1/4).*pi.^(-1).*((-1).*h.^(-1).*pi.*((-1).*cos(h.^(-1).*pi.*(y+...
            (-1).*y0))+cosh(h.^(-1).*pi.*(x+(-1).*x0))).^(-1).*sinh(h.^(-1).* ...
            pi.*(x+(-1).*x0))+(-1).*h.^(-1).*pi.*((-1).*cos(h.^(-1).*pi.*(y+( ...
            -1).*y0i))+cosh(h.^(-1).*pi.*(x+(-1).*x0))).^(-1).*sinh(h.^(-1).* ...
            pi.*(x+(-1).*x0)));
        
        Gx = Gx*mass_flow;
        
        Gy = (1/4).*pi.^(-1).*((-1).*h.^(-1).*pi.*((-1).*cos(h.^(-1).*pi.*(y+( ...
            -1).*y0))+cosh(h.^(-1).*pi.*(x+(-1).*x0))).^(-1).*sin(h.^(-1).* ...
            pi.*(y+(-1).*y0))+(-1).*h.^(-1).*pi.*((-1).*cos(h.^(-1).*pi.*(y+( ...
            -1).*y0i))+cosh(h.^(-1).*pi.*(x+(-1).*x0))).^(-1).*sin(h.^(-1).* ...
            pi.*(y+(-1).*y0i)));
        
        Gy = Gy*mass_flow;

%compute second derivative (first for velocity)

        %MATHEMATICA
        Gxx = (1/4).*pi.^(-1).*((-1).*h.^(-2).*pi.^2.*cosh(h.^(-1).*pi.*(x+(-1) ...
            .*x0)).*((-1).*cos(h.^(-1).*pi.*(y+(-1).*y0))+cosh(h.^(-1).*pi.*( ...
            x+(-1).*x0))).^(-1)+(-1).*h.^(-2).*pi.^2.*cosh(h.^(-1).*pi.*(x+( ...
            -1).*x0)).*((-1).*cos(h.^(-1).*pi.*(y+(-1).*y0i))+cosh(h.^(-1).* ...
            pi.*(x+(-1).*x0))).^(-1)+h.^(-2).*pi.^2.*((-1).*cos(h.^(-1).*pi.*( ...
            y+(-1).*y0))+cosh(h.^(-1).*pi.*(x+(-1).*x0))).^(-2).*sinh(h.^(-1) ...
            .*pi.*(x+(-1).*x0)).^2+h.^(-2).*pi.^2.*((-1).*cos(h.^(-1).*pi.*(y+ ...
            (-1).*y0i))+cosh(h.^(-1).*pi.*(x+(-1).*x0))).^(-2).*sinh(h.^(-1).* ...
            pi.*(x+(-1).*x0)).^2);
        
        Gxx = Gxx*mass_flow;
        
        Gxy = (1/4).*pi.^(-1).*(h.^(-2).*pi.^2.*((-1).*cos(h.^(-1).*pi.*(y+(-1) ...
            .*y0))+cosh(h.^(-1).*pi.*(x+(-1).*x0))).^(-2).*sin(h.^(-1).*pi.*( ...
            y+(-1).*y0)).*sinh(h.^(-1).*pi.*(x+(-1).*x0))+h.^(-2).*pi.^2.*(( ...
            -1).*cos(h.^(-1).*pi.*(y+(-1).*y0i))+cosh(h.^(-1).*pi.*(x+(-1).* ...
            x0))).^(-2).*sin(h.^(-1).*pi.*(y+(-1).*y0i)).*sinh(h.^(-1).*pi.*( ...
            x+(-1).*x0)));
        
        Gxy = Gxy*mass_flow;
        
        Gyx = (1/4).*pi.^(-1).*(h.^(-2).*pi.^2.*((-1).*cos(h.^(-1).*pi.*(y+(-1) ...
            .*y0))+cosh(h.^(-1).*pi.*(x+(-1).*x0))).^(-2).*sin(h.^(-1).*pi.*( ...
            y+(-1).*y0)).*sinh(h.^(-1).*pi.*(x+(-1).*x0))+h.^(-2).*pi.^2.*(( ...
            -1).*cos(h.^(-1).*pi.*(y+(-1).*y0i))+cosh(h.^(-1).*pi.*(x+(-1).* ...
            x0))).^(-2).*sin(h.^(-1).*pi.*(y+(-1).*y0i)).*sinh(h.^(-1).*pi.*( ...
            x+(-1).*x0)));
        
        Gyx = Gyx*mass_flow;
        
        Gyy = (1/4).*pi.^(-1).*((-1).*h.^(-2).*pi.^2.*cos(h.^(-1).*pi.*(y+(-1).* ...
            y0)).*((-1).*cos(h.^(-1).*pi.*(y+(-1).*y0))+cosh(h.^(-1).*pi.*(x+( ...
            -1).*x0))).^(-1)+(-1).*h.^(-2).*pi.^2.*cos(h.^(-1).*pi.*(y+(-1).* ...
            y0i)).*((-1).*cos(h.^(-1).*pi.*(y+(-1).*y0i))+cosh(h.^(-1).*pi.*( ...
            x+(-1).*x0))).^(-1)+h.^(-2).*pi.^2.*((-1).*cos(h.^(-1).*pi.*(y+( ...
            -1).*y0))+cosh(h.^(-1).*pi.*(x+(-1).*x0))).^(-2).*sin(h.^(-1).* ...
            pi.*(y+(-1).*y0)).^2+h.^(-2).*pi.^2.*((-1).*cos(h.^(-1).*pi.*(y+( ...
            -1).*y0i))+cosh(h.^(-1).*pi.*(x+(-1).*x0))).^(-2).*sin(h.^(-1).* ...
            pi.*(y+(-1).*y0i)).^2);
        
        Gyy = Gyy*mass_flow;
        
        %soultion of Stokes isung a potential field has constant prssure
        p = p_out;
        
        %Stress Tensor
        NXX = -p + 2*visc*Gxx;
        NXY = visc*Gxy + visc*Gyx;
        NYX = NXY;
        NYY = -p + 2*visc*Gyy;

%       cf2 = -pi/4/h^2;
% 
%       NUMx = cosh(A).*C-sinh(A).^2;
%       NUMxi = cosh(Ai).*Ci-sinh(Ai).^2;
%       Gxx = cf2 * (NUMx./C.^2 + NUMxi./Ci.^2);
%       
%       Gxy = 1;
%       
%       Gyx = Gxy;
%       
%       NUMy = cos(B).*C-sin(B).^2;
%       NUMyi = cos(Bi).*Ci-sin(Bi).^2;
%       Gyy = cf2 * (NUMy./C.^2 + NUMyi./Ci.^2);

%     figure
%     contour(x,y,G,100)
%     axis equal
%     xlabel('x')
%     ylabel('y')
%     title('velocity potential')
%     
%     figure
%     quiver(x,y,Gx,Gy)
%     axis equal
%     title('velocity field')
%     xlabel('x')
%     ylabel('y')
%     
%     figure
%     contour(x,y,p,100)
%     axis equal
%     xlabel('x')
%     ylabel('y')
%     title('pressure')

end
