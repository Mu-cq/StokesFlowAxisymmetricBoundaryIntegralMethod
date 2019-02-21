%remesh panels with distribution

function [xNew,yNew] = remeshDistrPanels2(t,x,y,PARAM,panel)

%error('Use remeshDistrPanels')

%compute arc lenght
dx = diff(x);   dy = diff(y);
dl = sqrt(dx.^2+dy.^2);
l = [0 cumsum(dl)];
l0 = sum(dl);
dlMax = max(dl);    dlMin = min(dl);
%RcurvMax = 1e3;

%kind of distribution
if PARAM.distr(panel)==0
    
    %uniform distribution
    distrHere = @(nElem) linspace(0,l0,nElem);
    
elseif PARAM.distr(panel)==1
    %cluster point where the curvature is higher
    
    %compute curvature
    [k1,k2] = computeCurvatureSplines(x,y,PARAM.orderVariableStokes(panel));
    %k = k1+k2;
    %k = k1;
    k = k1.*k1 + k2.*k2 + 0.004;
    kMesh = abs((k(1:end-1)+k(2:end))/2)+1;
    %RcurvMax = max(1./sqrt(k));
    
    %ditribution
    lMesh = ([0 cumsum((1./kMesh).^PARAM.adaptDistr(panel))]);  lMesh = lMesh/lMesh(end)*l0;
    pp = spline(linspace(0,l0,numel(lMesh)),lMesh);
    distrHere = @(nElem) (ppval(pp,linspace(0,l0,nElem)));
    
elseif PARAM.distr(panel)==2
    
    %cluster point cose to upper wall
    yMesh = (y(1:end-1)+y(2:end))/2;
    %distribution = (1./exp(yMesh)).^PARAM.adaptDistr(panel);
    distribution = exp(1-yMesh).^PARAM.adaptDistr(panel);
    lMesh = ([0 cumsum(distribution)]);  lMesh = lMesh/lMesh(end)*l0;
    pp = spline(linspace(0,1,numel(lMesh)),lMesh);
    distrHere = @(nElem) (ppval(pp,linspace(0,1,nElem)));
    
elseif PARAM.distr(panel)==3
    
    %cluster point cose to lateral wall
    xMesh = (x(1:end-1)+x(2:end))/2;
    distribution = exp(xMesh).^PARAM.adaptDistr(panel);
    lMesh = ([0 cumsum(distribution)]);  lMesh = lMesh/lMesh(end)*l0;
    pp = spline(linspace(0,1,numel(lMesh)),lMesh);
    distrHere = @(nElem) (ppval(pp,linspace(0,1,nElem)));
    
else
    error('Not implemented')
end

%check whether the distribution is respected
weight = integrationOnLineWeightAxis(x,y,PARAM.orderVariableStokes(panel),PARAM.orderGeometryStokes(panel),PARAM.SPlinesType(panel));
normCheck = weight*((l-distrHere(numel(l))).^2)';

xNew = x;   yNew = y;
if dlMax>PARAM.maxElem(panel) || dlMin<PARAM.minSizeElemRemesh(panel) || t==0 || normCheck>PARAM.normRemesh(panel)
    
    %apply distribution
    l = l/l(end)*l0;
    xNew = spline(l,x,distrHere(numel(x)));
    yNew = spline(l,y,distrHere(numel(x)));
    dx = diff(xNew);   dy = diff(yNew);
    dl = sqrt(dx.^2+dy.^2);
    l = [0 cumsum(dl)];
    l = l/l(end)*l0;
    dlMax = max(dl);    dlMin = min(dl);

    if dlMax>PARAM.maxElem(panel)

        %increase number of elements
        while dlMax>PARAM.maxElem(panel)
                        
            %new number of elements
            nElem = ceil(numel(xNew)+0.01*numel(xNew));

            display(['Add ' num2str(nElem-numel(xNew)) ' nodes to bubble, elem=' num2str(nElem)])

            %add elements
            xNew = spline(l,xNew,distrHere(nElem));
            yNew = spline(l,yNew,distrHere(nElem));

            dx = diff(xNew);   dy = diff(yNew);
            dl = sqrt(dx.^2+dy.^2);
            l = [0 cumsum(dl)];
            l = l/l(end)*l0;
            dlMax = max(sqrt(diff(xNew).^2+diff(yNew).^2));
            dlMin = min(sqrt(diff(xNew).^2+diff(yNew).^2));
            
            if nElem>PARAM.maxNumberTotalElem
                disp('Too many elements')
                break;
            end

        end

    elseif dlMin<PARAM.minSizeElemRemesh(panel)

        %decrease number of elements
        while dlMin<PARAM.minSizeElemRemesh(panel)

            %new number of elements
            nElem = floor(numel(xNew)-0.01*numel(xNew));

            display(['Eliminate ' num2str(-nElem+numel(xNew)) ' nodes from bubble, elem=' num2str(nElem)])

            %add elements
            xNew = spline(l,xNew,distrHere(nElem));
            yNew = spline(l,yNew,distrHere(nElem));

            dx = diff(xNew);   dy = diff(yNew);
            dl = sqrt(dx.^2+dy.^2);
            l = [0 cumsum(dl)];
            l = l/l(end)*l0;
            dlMin = min(sqrt(diff(xNew).^2+diff(yNew).^2));

        end
        
    else
        
        disp('Redistribuite nodes')

    end

end