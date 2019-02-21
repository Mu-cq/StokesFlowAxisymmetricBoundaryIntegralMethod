%remesh using a chosen distribution, this also can change the number of
%elements: new number of elements is equal to the numbear od data points in
%the distribution

function [x, y] = remeshThickWall(distr,PARAM)

    %geometry properties
    l = PARAM.L;            % wall lenght
    R = PARAM.R;            % motor radius
    r = PARAM.thickness;    % radius of curved part
    theta = PARAM.theta;    % wall inclination
    
    %select lower straight part

    %length of every element
    %ds = diff(curvCoord);
    s = l+pi*r;                %length of the all line
    %curv = [0 cumsum(ds)];      %curvilinear coordinate
    
    int = cumsum(distr);
    space = int/int(end)*s;              %spacing for the nodes with new distribution
    %space = space/space(end)*s;     %spacing for the nodes with new distribution
     
    %decide for the spacing between the nodes (based on distribution)
    s = [0 space'];                                         % modified curvilinear coordinate
    x1 = 0;    y1 = 0;                                     % center of the small circle
    phi = (s-l)/r;                                 % angle for small circle
    x = (l-s).*(s<l) + (r.*cos(-phi+3*pi/2)+x1).*(s>=l);
    y = -r.*(s<l) + (r.*sin(-phi+3*pi/2)+y1).*(s>=l);
    
    %rotation
    x = (x-l/2)*cos(theta) - y*sin(theta) + l/2;
    y = (x-l/2)*sin(theta) + y*cos(theta) + R;

end