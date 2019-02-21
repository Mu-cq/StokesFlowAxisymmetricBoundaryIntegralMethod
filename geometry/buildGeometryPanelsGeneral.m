

function [x,y,PARAM,tParametricBase] = buildGeometryPanelsGeneral(PARAM)

%build parametric grid
tParametric = parametricGrid(PARAM);

%build physical shape
[x,y] = buildGeometryPanelsParametric(tParametric,PARAM);

%rigid body rotation
for i = 1:numel(PARAM.n)
    
    if PARAM.rotate(i)~=0
        if isnan(PARAM.xCrotate(i)) || isnan(PARAM.yCrotate(i))
            error('Provide pivoting point')
        end
        [x,y,tParametric,PARAM.x0_Circle,PARAM.y0_Circle,PARAM.xStart,PARAM.xEnd,PARAM.yStart,PARAM.yEnd,PARAM.thetaStart,PARAM.thetaEnd] = rigidBodyRotation(x,y,tParametric,PARAM.xCrotate(i),PARAM.yCrotate(i),PARAM.rotate(i),PARAM.x0_Circle,PARAM.y0_Circle,PARAM.xStart,PARAM.xEnd,PARAM.yStart,PARAM.yEnd,PARAM.thetaStart,PARAM.thetaEnd,PARAM.geometryPanel,i);
    end
    
end
tParametricBase = tParametric;