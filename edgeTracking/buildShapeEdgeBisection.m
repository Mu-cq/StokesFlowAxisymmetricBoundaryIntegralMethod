%compute new shape by bisection

function [Tstart,initialShape] = buildShapeEdgeBisection(Tedge,Yedge,PARAM,V0)
    
%check which is last shape in and out the bassin
[T1,Y1,T2,Y2] = checkInOrOutTheBassin(Tedge,Yedge,PARAM);

%average shapes at the right time, applying volume consevation
[Tstart,initialShape] = getShapesAndAverage(T1,Y1,T2,Y2,PARAM,V0);

if PARAM.remeshStart==1
    %execute remesh spectrally
    disp('Remesh spectrally')
    initialShape = remeshDropSpectral(initialShape,PARAM);
end