%split elements when one elemnt is close to the other

function tParametric = remeshSplitProximity2(x,y,tParametric,tParametricBase,PARAM,panel1,panel2)

%compute physical panel lenght
lPanel = panelLenght(x,y,tParametric,panel1,PARAM);

%compute minimum distance for every elemnt
dl = diff(lPanel);
x1 = x{panel1};     y1 = y{panel1};
dist = distWallDrop(x{panel2},y{panel2},(x1(1:end-1)+x1(2:end))/2,(y1(1:end-1)+y1(2:end))/2);
checkActiveRemesh = sum(dl>dist/PARAM.coeffDist);

%minimum of minimum
distance = min(dist);

%check there is need to remesh
if checkActiveRemesh>=1 && distance<PARAM.distActivateRemesh(panel1)
    
    %split element if needed
    tParametric{panel1} = splitElements(tParametricBase,panel1,panel2,PARAM);
    
elseif distance>=PARAM.distActivateRemesh(panel1)
    
    %when the two panels are far apart, use first mesh
    tParametric{panel1} = tParametricBase{panel1};
    
end
        
















