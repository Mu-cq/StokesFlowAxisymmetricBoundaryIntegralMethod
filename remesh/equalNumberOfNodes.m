

function [x1,y1,x2,y2] = equalNumberOfNodes(x1,y1,x2,y2,panel,PARAM)
%equal number of nodes of two different shapes

if numel(x1)>numel(x2)
    
    x = x2;
    y = x2;
    nNew = numel(x1);
    
elseif numel(x1)<numel(x2)
    
    x = x1;
    y = y1;
    nNew = numel(x2);
    
else
    
    error('Not implemented')
    
end

%compute arc lenght
dx = diff(x);   dy = diff(y);
dl = sqrt(dx.^2+dy.^2);
l = [0 cumsum(dl)];
l0 = sum(dl);

%kind of distribution
if PARAM.distr(panel)==0
    
    %uniform distribution
    distrHere = @(nElem) linspace(0,l0,nElem);
    
elseif PARAM.distr(panel)==1
    %cluster point where the curvature is higher
    
    %compute curvature
    [k1,k2] = computeCurvatureSplines(x,y,PARAM.orderVariableStokes(panel));
    k = k1+k2;
    kMesh = abs((k(1:end-1)+k(2:end))/2)+1;
    
    %distrHere = @(nElem,lEnd) linspace(0,lEnd,nElem);
    lMesh = ([0 cumsum((1./kMesh).^PARAM.adaptDistr(panel))]);  lMesh = lMesh/lMesh(end)*l0;
    pp = spline(linspace(0,l0,numel(lMesh)),lMesh);
    distrHere = @(nElem) (ppval(pp,linspace(0,l0,nElem)));
    
    %error('Validate more carefully')
    
elseif PARAM.distr(panel)==2
    
    %cluster point cose to upper wall
    yMesh = (y(1:end-1)+y(2:end))/2;
    %distribution = (1./exp(yMesh)).^PARAM.adaptDistr(panel);
    distribution = exp(1-yMesh).^PARAM.adaptDistr(panel);
    lMesh = ([0 cumsum(distribution)]);  lMesh = lMesh/lMesh(end)*l0;
    pp = spline(linspace(0,1,numel(lMesh)),lMesh);
    distrHere = @(nElem) (ppval(pp,linspace(0,1,nElem)));
    
    minDist = @(xNew,yNew) min(1-yNew) > min(sqrt(diff(xNew).^2+diff(yNew).^2));
    
    %pp = spline(l/l(end),lMesh);
    %distrHere = @(nElem) (ppval(pp,linspace(0,1,nElem)));
    
elseif PARAM.distr(panel)==3
    
    %cluster point cose to lateral wall
    xMesh = (x(1:end-1)+x(2:end))/2;
    distribution = exp(xMesh).^PARAM.adaptDistr(panel);
    lMesh = ([0 cumsum(distribution)]);  lMesh = lMesh/lMesh(end)*l0;
    pp = spline(linspace(0,1,numel(lMesh)),lMesh);
    distrHere = @(nElem) (ppval(pp,linspace(0,1,nElem)));
    
    minDist = @(xNew,yNew) min(xNew) > min(sqrt(diff(xNew).^2+diff(yNew).^2));
    
    %pp = spline(l/l(end),lMesh);
    %distrHere = @(nElem) (ppval(pp,linspace(0,1,nElem)));
    
else
    error('Not implemented')
end

%new shape
xNew = spline(l,x,distrHere(nNew));
yNew = spline(l,y,distrHere(nNew));

if numel(x1)>numel(x2)
    
    x2 = xNew;
    y2 = yNew;
    
elseif numel(x1)<numel(x2)
    
    x1 = xNew;
    y1 = yNew;
    
end




