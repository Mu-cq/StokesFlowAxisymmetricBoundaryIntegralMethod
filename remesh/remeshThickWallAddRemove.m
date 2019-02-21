%remesh using a chosen distribution, this also can change the number of
%elements: new number of elements is equal to the numbear od data points in
%the distribution

function [x, y, add, remove] = remeshThickWallAddRemove(x,y,dist,minimum,PARAM)

    %min element lenght
    minElem = PARAM.minElemWall;

    %critical distance is the sistance to the droplet interface
    critical = dist;

    %geometry properties
    l = PARAM.L;            % wall lenght
    R = PARAM.R;            % motor radius
    %r = PARAM.thickness;    % radius of curved part
    theta = PARAM.theta;    % wall inclination
    
    %center and unturn
%     x = x-(max(x)+min(x))/2;
%     y = y-R;
%     x = x*cos(-theta) - y*sin(-theta);
%     y = x*sin(-theta) + y*cos(-theta);
    
    %check if distance is less than critical
    add = 0;
    remove = 0;
    dl = sqrt(diff(x).^2+diff(y).^2);
    count = 1;
    for i = 1:numel(x)-1
        
        %add point
        if dl(i)>critical(i) && dl(i)>2*minElem
            
            x = [x(1:count) (x(count)+x(count+1))/2 x(count+1:end)];
            y = [y(1:count) (y(count)+y(count+1))/2 y(count+1:end)];
            
            add = add+1;
            count = count+1;
            
        %remove only if it was not original
        elseif critical(i)>PARAM.R/4 && dl(i)<minimum %&& (1-sum(x(count)==PARAM.xWallStart))
            
            x = [x(1:count-1) x(count+1:end)];
            y = [y(1:count-1) y(count+1:end)];
            
            remove = remove+1;
            count = count-1;
            
        end
        count = count+1;
        
    end
    
    %check if element is too long compared to the one after
    dl = sqrt(diff(x).^2+diff(y).^2);
    count = 1;
    for i = 1:numel(x)-2
        
        if dl(i)>2*dl(i+1)+0.1*dl(i+1)
            
            x = [x(1:count) (x(count)+x(count+1))/2 x(count+1:end)];
            y = [y(1:count) (y(count)+y(count+1))/2 y(count+1:end)];
            
            add = add+1;
            count = count+1;
            
        end
        count = count+1;
        
    end
    
    %check if element is too long compared to the one before
    dl = sqrt(diff(x).^2+diff(y).^2);
    count = 2;
    for i = 2:numel(x)-1
        
        if dl(i)>2*dl(i-1)+0.1*dl(i-1)
            
            x = [x(1:count) (x(count)+x(count+1))/2 x(count+1:end)];
            y = [y(1:count) (y(count)+y(count+1))/2 y(count+1:end)];
            
            add = add+1;
            count = count+1;
            
        end
        count = count+1;
        
    end
    
    %rotation
%     x = (x-l/2)*cos(theta) - y*sin(theta) + l;
%     y = (x-l/2)*sin(theta) + y*cos(theta) + R;

end