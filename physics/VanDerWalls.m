%compute Van der Walls forces from 2 to 1

function [fx,fy] = VanDerWalls(x1,y1,x2,y2,dl2,elem)

%     fx = zeros(numel(x1),1);
%     fy = zeros(numel(x1),1);
%     for i = 1:numel(x1)
%         
%         x1now = x1(i);
%         y1now = y1(i);
%         for k = 1:numel(x2)
% 
%             x2now = x2(k);
%             y2now = y2(k);
%             
%             %compute distance between the 2 nodes
%             dist = sqrt((x1now-x2now)^2+(y1now-y2now)^2);
%             %force direction
%             nx = (x1now-x2now)/sqrt((x1now-x2now)^2+(y1now-y2now)^2); ny = (y1now-y2now)/sqrt((x1now-x2now)^2+(y1now-y2now)^2);
%             
%             fx(i) = fx(i) + dist^(-3)*nx;
%             fy(i) = fy(i) + dist^(-3)*ny;
%         
%         end
%     end

    %dl2 = sqrt(diff(x2).^2+diff(y2).^2);
    
    fx = zeros(numel(x1),1);
    fy = zeros(numel(x1),1);
    for i = 1:numel(x1)
        
        x1now = x1(i);
        y1now = y1(i);
        
        %compute distance between the 2 nodes
        dist = sqrt((x1now-x2).^2+(y1now-y2).^2);
        %force direction
        nx = (x1now-x2)./sqrt((x1now-x2).^2+(y1now-y2).^2); ny = (y1now-y2)./sqrt((x1now-x2).^2+(y1now-y2).^2);
            
        %decide intergation constan or linear element
        if elem == 2
            fx(i) = 2*pi*sum((dist(1:end-1).^(-3).*nx(1:end-1).*y2(1:end-1)+dist(2:end).^(-3).*nx(2:end).*y2(2:end)).*dl2)/2;
            fy(i) = 2*pi*sum((dist(1:end-1).^(-3).*ny(1:end-1).*y2(1:end-1)+dist(2:end).^(-3).*ny(2:end).*y2(2:end)).*dl2)/2;
        elseif elem == 1
            fx(i) = 2*pi*sum(dist.^(-3).*nx.*dl2.*y2);
            fy(i) = 2*pi*sum(dist.^(-3).*ny.*dl2.*y2);
        end
        
%         for k = 1:numel(x2)
% 
%             x2now = x2(k);
%             y2now = y2(k);
%             
%             %compute distance between the 2 nodes
%             dist = sqrt((x1now-x2now)^2+(y1now-y2now)^2);
%             %force direction
%             nx = (x1now-x2now)/sqrt((x1now-x2now)^2+(y1now-y2now)^2); ny = (y1now-y2now)/sqrt((x1now-x2now)^2+(y1now-y2now)^2);
%             
%             fx(i) = fx(i) + dist^(-3)*nx;
%             fy(i) = fy(i) + dist^(-3)*ny;
%         
%         end
    end

end