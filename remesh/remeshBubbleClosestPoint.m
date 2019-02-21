%remesh panels with distribution

function [xNew,yNew] = remeshBubbleClosestPoint(t,x,y,PARAM,panel)

%kind of distribution
if PARAM.distr(panel)==0
    %uniform distribution
    
    error('Not implemented')
elseif PARAM.distr(panel)==1
    %cluster point where the curvature is higher
    
    error('Not implemented')
    
elseif PARAM.distr(panel)==2
    
    %cluster point cose to upper wall
    yMesh = (y(1:end-1)+y(2:end))/2;
    lMesh = ([0 cumsum((1./yMesh).^PARAM.adaptDistr(panel))]);  lMesh = lMesh/lMesh(end);
    pp = spline(linspace(0,1,numel(lMesh)),lMesh);
    distrHere = @(nElem,lEnd) (ppval(pp,linspace(0,1,nElem))*lEnd);
    
    %find closest point to wall
    
    
else
    error('Not implemented')
end

%compute arc lenght
dx = diff(x);   dy = diff(y);
dl = sqrt(dx.^2+dy.^2);
l = [0 cumsum(dl)];
dlMax = max(dl);    dlMin = min(dl);

xNew = x;   yNew = y;
if dlMax>PARAM.maxElem(panel) || dlMin<PARAM.minSizeElemRemesh(panel) || t==0
    
    %apply distribution
    xNew = spline(l,x,distrHere(numel(x),l(end)));
    yNew = spline(l,y,distrHere(numel(x),l(end)));
    dx = diff(xNew);   dy = diff(yNew);
    dl = sqrt(dx.^2+dy.^2);
    l = [0 cumsum(dl)];
    dlMax = max(dl);    dlMin = min(dl);

    if dlMax>PARAM.maxElem(panel)

        %increase number of elements
        while dlMax>PARAM.maxElem(panel)

            %new number of elements
            nElem = ceil(numel(xNew)+0.01*numel(xNew));

            display(['Add ' num2str(nElem-numel(xNew)) ' nodes to bubble, elem=' num2str(nElem)])

            %add elements
            xNew = spline(l,xNew,distrHere(nElem,l(end)));
            yNew = spline(l,yNew,distrHere(nElem,l(end)));

            dx = diff(xNew);   dy = diff(yNew);
            dl = sqrt(dx.^2+dy.^2);
            l = [0 cumsum(dl)];
            dlMax = max(sqrt(diff(xNew).^2+diff(yNew).^2));

        end

    elseif dlMin<PARAM.minSizeElemRemesh(panel)

        %decrease number of elements
        while dlMin<PARAM.minSizeElemRemesh(panel)

            %new number of elements
            nElem = floor(numel(xNew)-0.01*numel(xNew));

            display(['Eliminate ' num2str(-nElem+numel(xNew)) ' nodes from bubble, elem=' num2str(nElem)])

            %add elements
            xNew = spline(l,xNew,distrHere(nElem,l(end)));
            yNew = spline(l,yNew,distrHere(nElem,l(end)));

            dx = diff(xNew);   dy = diff(yNew);
            dl = sqrt(dx.^2+dy.^2);
            l = [0 cumsum(dl)];
            dlMin = min(sqrt(diff(xNew).^2+diff(yNew).^2));

        end
        
    else
        
        display('Redistribuite nodes')

    end

end