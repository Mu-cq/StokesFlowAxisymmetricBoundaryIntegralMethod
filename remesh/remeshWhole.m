%remesh first part of panel

function tParametric = remeshWhole(tParametric,panel,PARAM,distELement,distance)

%adaptivity coeff
coeff = PARAM.adaptCoeff(panel);
numelNodes = numel(tParametric);

tParametricOld = tParametric-tParametric(1);

%compute ditribution of dt size
distr = 1./distELement.^coeff;
%cumulative = [0 cumsum(distr)]/sum(distr);
dt = diff(tParametric);
cumulative = [0 cumsum(distr(1:end-1).*dt)]/sum(distr(1:end-1).*dt);

%remesh first half
count1 = 0;
minT = 1;
maxT = 1;
while maxT>PARAM.maxElem(panel) || minT>distance
        
        count1 = count1+1;
    
        if count1>1
            numelNodes = round(1.1*numelNodes);
        end
        
        tParametric = spline(cumulative,tParametricOld,linspace(tParametricOld(1),tParametricOld(end),numelNodes-1));
        
        if PARAM.geometryPanel(panel)==1 % is and arc
            tParametric = tParametric*(PARAM.thetaEnd(panel)-PARAM.thetaStart(panel))+PARAM.thetaStart(panel);
        end
        
        maxT = max(diff(tParametric))*2;
        minT = min(diff(tParametric));

        %physical space
        if PARAM.geometryPanel(panel)==1 % is and arc

            maxT = PARAM.rArc(panel)*maxT;
            minT = PARAM.rArc(panel)*minT;

        elseif PARAM.geometryPanel(panel)==0 % is a line

            dx = PARAM.xEnd(panel)-PARAM.xStart(panel);
            dy = PARAM.yEnd(panel)-PARAM.yStart(panel);
            dl = sqrt(dx^2+dy^2);
            maxT = dl*maxT;
            minT = dl*minT;

        else

            error('Not implemented')

        end
        
        if minT<2*PARAM.minSizeElemRemesh(panel)
            warning('Stop remesh because smallest elment is becoming too small')
            break
        end
        
        if count1>100
            warning('Too many iterations in remeshing')
            break
        end
    
end