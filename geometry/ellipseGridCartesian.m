%compute cartesina coordinates of an ellipse of given D on given grid

function [x,y] = ellipseGridCartesian(theta,D)

    %compute ellipse having volume as a unit raiuds sphere
    bAxis = ((1.0-D)/(1.0+D))^(1/3);
    aAxis = 1.0/bAxis^2;
    
    %compute cartesian coordinates
    x = aAxis*cos(theta);    y = bAxis*sin(theta);

end