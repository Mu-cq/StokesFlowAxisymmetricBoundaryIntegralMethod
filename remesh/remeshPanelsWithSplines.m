%remesh block following different algorithms

function varNew = remeshPanelsWithSplines(t,var,PARAM)

if t>0

%number of panels
numPanels = numel(PARAM.n);

%loop on panels
for i = 1:numPanels
    
   if PARAM.remeshType(i)==4
       
       %this is ook for one bubble with moving wall
       xDrop = var(1:2:end-2);
       yDrop = var(2:2:end-1);
       
       %remesh with uniform distribution, taking care that elements are not
       %too large or too small
       [xDrop,yDrop] = remeshDistrPanels(xDrop,yDrop,PARAM,i);
       var(1:2:end-2) = xDrop;
       var(2:2:end-1) = yDrop;
       varNew = var;
       
   end
   
end

end