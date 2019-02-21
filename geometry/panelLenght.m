%compute physical panel lenght

function lPanel = panelLenght(x,y,tParametric,panel1,PARAM)

if PARAM.geometryPanel(panel1)==0 % is a straight line
    
    xLine = x{panel1};
    yLine = y{panel1};
    
    dx = diff(xLine);   dy = diff(yLine);
    dl = sqrt(dx.^2+dy.^2);
    lPanel = [0 cumsum(dl)];
    
elseif PARAM.geometryPanel(panel1)==1 % is an arc
    
    rCircle = PARAM.rArc(panel1);
    tCircle = tParametric{panel1};
    dt = [0 cumsum(diff(tCircle))];
    
    lPanel = rCircle*dt;
    
else
    
    error('Not implemented')
    
end
