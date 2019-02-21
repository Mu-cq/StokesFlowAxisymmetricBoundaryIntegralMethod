%remesh panels during time stepping

function [x,y,tParametric,PARAM] = remeshPanelsFunction(t,x,y,tParametric,PARAM,tParametricBase,nPerLenght)

if t>0% && sum(i==1:PARAM.remeshStep:loop)
            tParametricBase{5} = linspace(0,pi,round(PARAM.rArc(5)*pi*nPerLenght));
            tParametric = remeshPanels(x,y,tParametric,tParametricBase,PARAM);
            for k = 1:numel(PARAM.n)
                PARAM.n(k) = numel(tParametric{k})-1;
            end
            %build shape
            [x,y] = buildGeometryPanelsParametric(tParametric,PARAM);
end