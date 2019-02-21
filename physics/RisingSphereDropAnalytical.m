%analytical solution of a spherical droplet rising in a viscous fluid (from LOW REYNOLDS NUMBER HYDRODYNAMIC section 4.21)

function [u_x,u_y,p] = RisingSphereDropAnalytical(x,y,visc,a,U,capillary)

    %viscosity ratio is definite the opposite in the book
    sigma = 1/visc;

    %make grid
    THETA = atan(y./x).*(x>0) + (atan(y./x)+pi).*(x<0) + pi/2*(x==0);
    R = sqrt(x.^2+y.^2);
    
    %constants
    C = 0.5*U;
    B = 1.5*U*a*(1+2/3*sigma)/(1+sigma);
    D = 0.25*U*a^3/(1+sigma);
    E = 2.5*U/a^2*sigma/(1+sigma);
    G = -0.25*U*sigma/(1+sigma);
    
    %streamlines for inner and outer domain
    PhiOut = sin(THETA).^2.*(-0.5*B*R + C*R.^2 + D./R).*(R>=a);
    PhiIn = sin(THETA).^2.*(0.1*E*R.^4 + G*R.^2).*(R<a);
    
    %streamlines
    Phi = PhiOut + PhiIn;
    
    %velocities as derivatives of the streamlines
    %r direction
    u_rOut = 2*cos(THETA).*(-0.5*B./R + C + D./R.^3).*(R>=a);
    u_rIn = 2*cos(THETA).*(0.1*E*R.^2 + G).*(R<a);
    u_r = u_rOut + u_rIn;
    %theta direction
    u_thetaOut = -sin(THETA).*(-0.5*B./R + 2*C - D./R.^3).*(R>=a);
    u_thetaIn = -sin(THETA).*(0.4*E*R.^2 + 2*G).*(R<a);
    u_theta = u_thetaOut + u_thetaIn;
    
    %cartesian component, change sign to respect convention in the book
    u_x = u_r.*cos(THETA) - u_theta.*sin(THETA);
    u_y = - u_r.*sin(THETA) - u_theta.*cos(THETA);
    
    %pressure jump at the interface
    %Dp = 2/a/capillary + a*cos(THETA)*3*(1+1.5*visc)/(1+visc);
    Dp = 2/capillary/a;
    %Dp = 2/capillary/a + a*cos(THETA);
    
    %p from section 4-15
    p_out = -B*cos(THETA)./R.^2.*(R>=a);
    %p_out = -2*B*cos(THETA)./R.*(R>=a);
    p_in = (-2*B*cos(THETA)./a^2 + Dp - 2*visc*E*cos(THETA).*(R-a)).*(R<a);
    p = - p_in - p_out;
    
    %plot streamlines
%     figure
%     contour(y,x,Phi,100)
%     hold on
%     contour(-y,x,Phi,100)
%     plot(a*sin(0:0.001:2*pi),a*cos(0:0.001:2*pi),'k','LineWidth',2)
%     hold off
%     axis equal
%     axis([-max(max(y)) max(max(y)) min(min(x)) max(max(x))])
%     grid on
%     xlabel('radial')
%     ylabel('axial')
%     title(['Streamlines \lambda=' num2str(visc) ', U=' num2str(U)])
%     
%     %plot velocity field
%     figure
%     quiver(y,x,u_y,u_x)
%     hold on
%     quiver(-y,x,-u_y,u_x,'b')
%     plot(a*sin(0:0.001:2*pi),a*cos(0:0.001:2*pi),'k','LineWidth',2)
%     hold off
%     axis equal
%     axis([-max(max(y)) max(max(y)) min(min(x)) max(max(x))])
%     %grid on
%     xlabel('radial')
%     ylabel('axial')
%     title(['Velocity fiels \lambda=' num2str(visc) ', U=' num2str(U)])

end