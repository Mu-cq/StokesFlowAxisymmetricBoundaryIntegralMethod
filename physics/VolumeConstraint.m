%compute velocity normal to the interface in an extensional flow

function [DisEquality,Equality] = VolumeConstraint(theta,r,Vfinal)
    
    %trapezi integral
    INT = ([diff(theta),0]+[0,diff(theta)])/2;
    
    %volume constraint
    DisEquality = [];
    Equality = (2/3*pi*INT*(r.^3.*sin(theta'))-Vfinal)';

end