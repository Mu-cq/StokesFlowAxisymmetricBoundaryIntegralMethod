%remesh using a chosen distribution, this also can change the number of
%elements: new number of elements is equal to the numbear od data points in
%the distribution

function [x, y, add, remove] = remeshDropAddRemove(x,y,dist,minimum,PARAM)

    %min element lenght
    minElem = PARAM.minElemDrop;

    %compute spline coeff
    [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(x, y);

    %critical distance is the sistance to the droplet interface
    critical = dist;
    
    %check if distance is less than critical
    add = 0;
    remove = 0;
    dl = sqrt(diff(x).^2+diff(y).^2);
    count = 1;
    for i = 1:numel(x)-1
        
        %add point
        if (dl(i)>critical(i) || dl(i)>2.5*minimum) && dl(i)>2*minElem
            
            %compute midlle point with splines
            t = 0.5;
            xMiddle = ax(count) + bx(count)*t + cx(count)*t^2 + dx(count)*t^3;
            yMiddle = ay(count) + by(count)*t + cy(count)*t^2 + dy(count)*t^3;
            
            x = [x(1:count) xMiddle x(count+1:end)];
            y = [y(1:count) yMiddle y(count+1:end)];
            
            [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(x, y);
            
            add = add+1;
            count = count+1;
            
        %remove
        elseif (critical(i)>PARAM.R/4 && dl(i)<minimum) || (dl(i)<0.5*critical(i) && dl(i)<minimum) %&& (1-sum(x(count)==PARAM.xWallStart))
            
            x = [x(1:count-1) x(count+1:end)];
            y = [y(1:count-1) y(count+1:end)];
            
            [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(x, y);
            
            remove = remove+1;
            count = count-1;
            
        end
        count = count+1;
        
    end
    
    %check if element is too long compared to the one after
%     dl = sqrt(diff(x).^2+diff(y).^2);
%     count = 1;
%     for i = 1:numel(x)-2
%         
%         if dl(i)>2*dl(i+1)+0.1*dl(i+1)
%             
%             %compute midlle point with splines
%             t = 0.5;
%             xMiddle = ax(count) + bx(count)*t + cx(count)*t^2 + dx(count)*t^3;
%             yMiddle = ay(count) + by(count)*t + cy(count)*t^2 + dy(count)*t^3;
%             
%             x = [x(1:count) xMiddle x(count+1:end)];
%             y = [y(1:count) yMiddle y(count+1:end)];
%             
%             [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(x, y);
%             
%             add = add+1;
%             count = count+1;
%             
%         end
%         count = count+1;
%         
%     end
%     
    %check if element is too long compared to the one before
    dl = sqrt(diff(x).^2+diff(y).^2);
    count = 2;
    for i = 2:numel(x)-1
        
        if dl(i)>2*dl(i-1)+0.1*dl(i-1)
            
            %compute midlle point with splines
            t = 0.5;
            xMiddle = ax(count) + bx(count)*t + cx(count)*t^2 + dx(count)*t^3;
            yMiddle = ay(count) + by(count)*t + cy(count)*t^2 + dy(count)*t^3;
            
            x = [x(1:count) xMiddle x(count+1:end)];
            y = [y(1:count) yMiddle y(count+1:end)];
            
            [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(x, y);
            
            add = add+1;
            count = count+1;
            
        end
        count = count+1;
        
    end
    
    display(['Add ' num2str(add) ' elements to the drop'])
    display(['Remove ' num2str(remove) ' elements from the drop'])

end