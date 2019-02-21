%remesh using a chosen distribution, this also can change the number of
%elements: new number of elements is equal to the numbear od data points in
%the distribution

function [x, y, add, remove] = remeshCircleLeftAddRemove(x,y,dist,minimum,PARAM)

    %critical distance is the sistance to the droplet interface
    critical = min(dist)*2;

    %geometry properties
    l = PARAM.L;            % wall lenght
    R = PARAM.R;            % motor radius
    r = PARAM.thickness;    % radius of curved part
    theta = PARAM.theta;    % wall inclination
    
    %center and unturn
    x = x-(max(x)+min(x))/2;
    y = y-R;
    x = x*cos(-theta) - y*sin(-theta);
    y = x*sin(-theta) + y*cos(-theta);
    
    %theta of point of this circle
    thetaHere = atan(y./x) + pi*(x<0).*(y>0) + pi*(x<0).*(y<0) + 2*pi*(x>0).*(y<0);
    
    %figure out if it left or right circle
    sign = (min(x)>0) - (min(x)<0);
    
    %check if distance is less than critical
    add = 0;
    remove = 0;
    dl = sqrt(diff(x).^2+diff(y).^2);
    for i = 1:numel(x)-1
        
        %add point
        if dist(i)<critical
            
            thetaNew = (thetaHere(i)+thetaHere(i+1))/2;
            x = [x(1:i) r*cos(thetaNew)+sign*L/2 x(i+1:end)];
            y = [y(1:i) r*sin(thetaNew) y(i+1:end)];
            
            add = add+1;
            
        %remove
        elseif dist(i)>4*critical&&dl(i)<minimum
            
            x = [x(1:i-1) x(i+1:end)];
            y = [y(1:i-1) y(i+1:end)];
            
            remove = remove+1;
            
        end
        
    end
    
    %rotation
    x = (x-l/2)*cos(theta) - y*sin(theta) + l/2;
    y = (x-l/2)*sin(theta) + y*cos(theta) + R;

end