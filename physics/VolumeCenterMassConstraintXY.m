%compute velocity normal to the interface in an extensional flow

function [DisEquality,Equality] = VolumeCenterMassConstraintXY(XY,Vfinal,xcmFinal)
    
    %n = numel(XY)/2;

    %cartesian coordiantes
    x = XY([1 2:2:end]);
    y = [0 XY(3:2:end-1) 0];

    %volume constraint
    DisEquality = [];
    EqualityVol = axis_int_gauss_vect(x,y)-Vfinal;
    EqualityCM = center_mass(x,y) -xcmFinal;
    
    Equality = [EqualityVol EqualityCM];

end