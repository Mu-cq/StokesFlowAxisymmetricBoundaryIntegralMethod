%compute surface integral with trapezi rule in cylindrical coordiante

function [A,dA] = surf_lailai_simple(x,y)

    %using trapezi rules for integration
    %dtheta = pi/(numel(x)-1);
    
    xcm = center_mass(x,y);
    x = x-xcm;
    
    R = sqrt(x.*x+y.*y);
        
    theta = acos(x./R);
    dtheta = theta(2:end)-theta(1:end-1);
    
%     hold on
%     plot(R.*cos(theta),R.*sin(theta),'or-')
%     hold off

%     figure
%     plot(R)
    
    sintheta = y./R;
    dRdtheta = [diff(R)./diff(theta) diff(R(end-1:end))./diff(theta(end-1:end))];
    
    tot = R.*sintheta.*sqrt(dRdtheta.*dRdtheta+R.*R);
    %tot = R.*R;
    dA = pi*(tot(1:end-1)+tot(2:end)).*dtheta;
    
    A = sum(dA);
    %A = trapz(theta,tot);   
    

end