%remesh block following different algoruthms

function tParametric = remeshPanels(x,y,tParametric,tParametricBase,PARAM)

%number of panels
numPanels = numel(PARAM.n);

%loop on panels
for i = 1:numPanels
    
   if PARAM.remeshType(i)==1 || PARAM.remeshType(i)==3
       
       %display('Remesh is activated, split elements into two')
       whichPanelIsClose = PARAM.remeshProximity{i};
       
       %find out which is the closest panel
       x1 = x{i};   y1 = y{i};
       minDist = zeros(numel(whichPanelIsClose),1);
       for k = 1:numel(whichPanelIsClose)
       
           x2 = x{whichPanelIsClose(k)};   y2 = y{whichPanelIsClose(k)};
           dist = distWallDrop(x1,y1,x2,y2);
           minDist(k) = min(dist);
       
       end
       
       [~,minInd] = min(minDist);
       if PARAM.remeshType(i)==1
           tParametric = remeshSplitProximity2(x,y,tParametric,tParametricBase,PARAM,i,whichPanelIsClose(minInd));
       elseif PARAM.remeshType(i)==3
           error('Bug')
           tParametric = remeshSplitProximity3(x,y,tParametric,PARAM,i,whichPanelIsClose(minInd));
       end
       
   elseif PARAM.remeshType(i)==2
       
       error('Bug')
       
       %remesh using distribution
       whichPanelIsClose = PARAM.remeshProximity{i};
       
       %find out which is the closest panel
       x1 = x{i};   y1 = y{i};
       minDist = zeros(numel(whichPanelIsClose),1);
       for k = 1:numel(whichPanelIsClose)
       
           x2 = x{k};   y2 = y{k};
           dist = distWallDrop(x1,y1,x2,y2);
           minDist(k) = min(dist);
       
       end
       
       [~,minInd] = min(minDist);
       tParametric = remeshDistribution(x,y,tParametric,PARAM,i,whichPanelIsClose(minInd));
    
   end
   
end