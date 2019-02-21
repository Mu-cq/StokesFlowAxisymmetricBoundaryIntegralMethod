%compute new shape by bisection

function [Tstart,initialShape] = buildShapeEdgeExtrapolateBEM(Tedge,Yedge,PARAM,V0)
    
%check which is last shape in and out the bassin
[T1,Y1,T2,Y2] = checkInOrOutTheBassinBEM(Tedge,Yedge,PARAM);

%average shapes at the right time, applying volume consevation
[Tstart,initialShape] = getShapesAndExtrapolateBEM(T1,Y1,T2,Y2,PARAM,V0);