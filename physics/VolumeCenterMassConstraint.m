%compute velocity normal to the interface in an extensional flow

function [DisEquality,Equality] = VolumeCenterMassConstraint(theta,r,Vfinal,xcmFinal)
    
    %trapezi integral
    INT = ([diff(theta),0]+[0,diff(theta)])/2;
    
    %volume constraint
    DisEquality = [];
    EqualityVol = 2/3*pi*INT*(r.^3.*sin(theta'))-Vfinal;
    EqualityCM = 3/4 * (INT*(r.^4.*sin(theta').*cos(theta')))/(INT*(r.^3.*sin(theta'))) -xcmFinal;
    
    Equality = [EqualityVol EqualityCM];

end