%compute block coordinates


function [xHere,yHere,nxHere,nyHere,panelRange,dlHere] = getBlockCoordinates(x,y,PARAM,block)

%number of panels in this block
panels = PARAM.panels(block);

%previous panels
if block==1
    startPan = 0;
else
    startPan = sum(PARAM.panels(1:block-1));
end

panelRange = startPan+1:startPan+panels;
xHere = [];
yHere = [];
nxHere = [];
nyHere = [];
dlHere = [];
for i = panelRange
    
   xNow = x{i};
   yNow = y{i};
   
   %compute normal vector
   [nxNow,nyNow] = computeNormalVector(xNow,yNow,PARAM.orderVariable(i),PARAM.orderGeometry(i),PARAM.SPlinesType(i));
   
   %arc lenght on this block
   dlNow = computeDL(xNow,yNow,PARAM.orderVariable(i),PARAM.orderGeometry(i),PARAM.SPlinesType(i));
   dlHere = [dlHere dlNow];
   
   if i>panelRange(1)
       xNow = xNow(2:end);
       yNow = yNow(2:end);
       if PARAM.orderVariable(panelRange(i))~=0
           nxNow = nxNow(2:end);
           nyNow = nyNow(2:end);
       end
   end
    
   xHere = [xHere xNow];
   yHere = [yHere yNow];
   
   nxHere = [nxHere nxNow];
   nyHere = [nyHere nyNow];
    
end
