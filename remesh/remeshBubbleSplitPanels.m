%remesh panels with distribution

function [xNew,yNew,checkActiveRemesh] = remeshBubbleSplitPanels(t,x,y,PARAM,panel1)

if numel(PARAM.remeshProximity{panel1})>1
    error('This remesh works only in proximity to one specific panel')
else
    panel2 = PARAM.remeshProximity{panel1};
end

%compute distance between points and panel
if PARAM.geometryPanel(panel2)==0
   m = (PARAM.yEnd(panel2)-PARAM.yStart(panel2))/(PARAM.xEnd(panel2)-PARAM.xStart(panel2));
   aaa = m;
   bbb = -1;
   ccc = PARAM.yStart(panel2)-m*PARAM.xStart(panel2);
   
   xm = (x(1:end-1)+x(2:end))/2;
   ym = (y(1:end-1)+y(2:end))/2;
   dist = abs(aaa*xm+bbb*ym+ccc)/sqrt(aaa^2+bbb^2);
else
   error('Not implemented')
end

%compute minimum distance for every element
dx = diff(x);   dy = diff(y);   dl = sqrt(dx.^2+dy.^2);
checkActiveRemesh = sum(dl>dist/PARAM.coeffDist);

%minimum of minimum
distance = min(dist);

%check there is need to remesh
if checkActiveRemesh>=1 && distance<PARAM.distActivateRemesh(panel1) && isempty(t)==0
    
    %split element if needed
    [xNew,yNew] = splitElementsBubble(x,y,panel1,dist,PARAM);
    
else
    
    xNew = x;
    yNew = y;
    
end





