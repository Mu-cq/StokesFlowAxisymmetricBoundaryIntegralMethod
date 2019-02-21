%replace interface point too close to walls

function [xDrop,yDrop] = replaceInterface(xDrop,yDrop,PARAM,distDrop)
    
    %compute spline coeff
    [~, bx, cx, dx, ~, by, cy, dy] = spline_symmetric(xDrop, yDrop);
    
    %compute normal vector
    N = [by./sqrt(bx.*bx+by.*by) (by(end)+2*cy(end)+3*dy(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)));...
      -bx./sqrt(bx.*bx+by.*by) (-bx(end)-2*cx(end)-3*dx(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)))];
    
%     for i = 1:numel(xDrop)
%         
%         if distDrop(i)<PARAM.repulsiveOn
%             
%             display('displace')
%            
%             %amount of displacment needed
%             displace = PARAM.repulsiveOn-distDrop(i);
%             
%             %diplace point along itts normal
%             xDrop(i) = xDrop(i) - displace*N(1,i);
%             yDrop(i) = yDrop(i) - displace*N(2,i);
%             
%         end
%         
%     end
    
    %display('displace')

    %amount of displacment needed
    displace = PARAM.repulsiveOn-distDrop;
    
    %diplace point along itts normal
    xDrop = xDrop - displace.*N(1,:).*(distDrop<PARAM.repulsiveOn);
    yDrop = yDrop - displace.*N(2,:).*(distDrop<PARAM.repulsiveOn);

end