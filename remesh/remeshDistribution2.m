%remesh analytical shape (panel 1) using distribution which might be based
%upon the proximity to panel 2

function tParametric = remeshDistribution2(x,y,tParametric,PARAM,panel1,panel2)

%compute distance between two panels
[distance,indNode] = panelDistance(x,y,panel1,panel2,PARAM);

%compute physical panel lenght
lPanel = panelLenght(x,y,tParametric,panel1,PARAM);
%minDL = min(diff(lPanel));
if indNode==1
    indNode = indNode+1;
elseif indNode==numel(lPanel)
    indNode = indNode-1;
end
minDLcloseToOtherPanel = lPanel(indNode)-lPanel(indNode-1);

%compute first and second half parts
tParametricHere = tParametric{panel1};

%check there is need to remesh
if minDLcloseToOtherPanel>distance/2 && distance<PARAM.distActivateRemesh(panel1)
    
    %distance per element
    dist = distWallDrop(x{panel2},y{panel2},x{panel1},y{panel1});
    
    %remesh first half
    tParametricHere = remeshWhole(tParametricHere,panel1,PARAM,dist,distance);
    
    %merge first and second half
    tParametric{panel1} = tParametricHere;
    
end

%check overlapping
if abs(diff(tParametric{panel1}))<1e-5
        error('Nodes overlapping')
end










