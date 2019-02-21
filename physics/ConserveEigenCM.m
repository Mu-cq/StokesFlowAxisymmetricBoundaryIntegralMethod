%conserve volume via an eigenfunction

function R = ConserveEigenCM(theta,alpha,eigen,r,XCM)

    %integration
    %INT = ([diff(theta),0]+[0,diff(theta)])/2;
    
    r = r+alpha*eigen;
    x = r'.*cos(theta);
    y = r'.*sin(theta);

    R = center_mass(x,y) - XCM;
    
end