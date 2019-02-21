%conserve volume via an eigenfunction

function R = ConserveEigenVolumeCM(theta,alpha,eigen,r,V0,XCM)

    %each variable
    alpha1 = alpha(1);
    alpha2 = alpha(2);
    
    %each eigenfunction
    eigen1 = eigen(:,1);
    eigen2 = eigen(:,2);

    %integration
    INT = ([diff(theta),0]+[0,diff(theta)])/2;
    
    %center of mass
    rMass = r + alpha1*eigen1 + alpha2*eigen2;
    x = rMass'.*cos(theta);
    y = rMass'.*sin(theta);

    %volume conservation
    %R = (2/3*pi*INT*((r + alpha1*eigen1+alpha2*eigen2).^3.*sin(theta')) - V0)^2 + (center_mass(x,y) - XCM)^2;
    
    R(1) = 2/3*pi*INT*((r + alpha1*eigen1 + alpha2*eigen2).^3.*sin(theta')) - V0;
    R(2) = center_mass(x,y) - XCM;
    
    
end