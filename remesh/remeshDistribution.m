%remesh analytical shape (panel 1) using distribution which might be based
%upon the proximity to panel 2

function tParametric = remeshDistribution(x,y,tParametric,PARAM,panel1,panel2)

%compute distance between two panels
[distance,indNode,distanceNode] = panelDistance(x,y,panel1,panel2,PARAM);
% hold on
% xHere = x{panel1};
% yHere = y{panel1};
% plot(xHere(indNode),yHere(indNode),'ok')

%compute physical panel lenght
lPanel = panelLenght(x,y,tParametric,panel1,PARAM);
%minDL = min(diff(lPanel));
if indNode==1
    indNode = indNode+1;
elseif indNode==numel(lPanel)
    indNode = indNode-1;
end
minDLcloseToOtherPanel = lPanel(indNode)-lPanel(indNode-1);
l0 = lPanel(indNode);
lEnd = lPanel(end);
if PARAM.geometryPanel(panel1)==1
    
    theta0 = PARAM.thetaStart(panel1);
    lenghtTheta = PARAM.thetaEnd(panel1)-PARAM.thetaStart(panel1);
    l0 = l0/lPanel(end)*lenghtTheta+theta0;
    lEnd = lEnd/lPanel(end)*lenghtTheta+theta0;
    
elseif PARAM.geometryPanel(panel1)==0
    
    lenghtL = lPanel(end);
    l0 = l0/lenghtL;
    lEnd = lEnd/lenghtL;
    
else
    
    error('Not implemented')
    
end

%compute first and second half parts
tParametricHere = tParametric{panel1};
tParametricFirst = tParametricHere(1:indNode);
tParametricSecond = tParametricHere(indNode:end);

%check there is need to remesh
if minDLcloseToOtherPanel>distance/2 && distance<PARAM.distActivateRemesh(panel1)
    
    %remesh first half
    tParametricFirst = remeshFirstPart(tParametricFirst,l0,panel1,PARAM,distance);
    
    %remesh second half
    tParametricSecond = remeshSecondPart(tParametricSecond,l0,lEnd,panel1,PARAM,distance);
    
    %merge first and second half
    tParametric{panel1} = [tParametricFirst(1:end-1) tParametricSecond];
    
% elseif distance>PARAM.distActivateRemesh(panel1)
%     
%     %error('test')
% 
%     %uniform distribution of given element size
%     nElem = round(lPanel(end)/PARAM.maxElem(panel1))+1;
%     tParametric{panel1} = linspace(0,1,nElem);
% 
%     if PARAM.geometryPanel(panel1)==1
% 
%         tParametric{panel1} = tParametric{panel1}*(PARAM.thetaEnd(panel1)-PARAM.thetaStart(panel1))+PARAM.thetaStart(panel1);
% 
%     end
    
end

%check overlapping
if abs(diff(tParametric{panel1}))<1e-5
        error('Nodes overlapping')
end










