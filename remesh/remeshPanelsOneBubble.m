%remesh block following different algorithms

function [varNew,status] = remeshPanelsOneBubble(t,var,PARAM)

status = 0;
if numel(t)==1

%number of panels
numPanels = numel(PARAM.n);

%loop on panels
for i = 1:numPanels
    
   %this is ok for one bubble with moving wall
   nNodes = floor(numel(var)/2);
   xDrop = var(1:2:2*nNodes-1)';
   yDrop = var(2:2:2*nNodes)';
    
   if PARAM.remeshType(i)==4
       
       %remesh with uniform distribution, taking care that elements are not
       %too large or too small
       [xDrop,yDrop] = remeshDistrPanels2(t,xDrop,yDrop,PARAM,i);
      
   elseif PARAM.remeshType(i)==5
       
       %remesh by splitting elements into two on spline
       [xDrop,yDrop] = remeshBubbleSplitPanels(t,xDrop,yDrop,PARAM,i);
       
   end
   
   varNew = zeros(2*numel(xDrop),1);
   varNew(1:2:2*numel(xDrop)-1) = xDrop;
   varNew(2:2:2*numel(xDrop)) = yDrop;
   
   if numel(xDrop)>PARAM.maxNumberTotalElem
       disp('Too many elements')
       status = 1;
   end
   
end

if mod(numel(var),2)==1
   varNew = [varNew; var(end)];
end

end