%give the solution in terms of stress and velocities
%having the following input: radius of the pipe R, length of the pipe L,
%number of element at the outlet n, number of element at the wall m, number of
%element at the inlet j, inlet velocity
%vel_in, vertical stress at the outlet "stress" and the viscosity of the
%fluid

function [xxx,Xsing,Ysing,nnx,nny,nnn] = BEM_Stokes(x,y,PARAM)

  %change variables name
  PARAM.typeBC = PARAM.typeBCstokes;
  PARAM.orderVariable = PARAM.orderVariableStokes;
  PARAM.orderGeometry = PARAM.orderGeometryStokes;
  
  %compute location of the singularity
  [Xsing,Ysing,nnn] = computeSingularityLocation(x,y,PARAM);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %compute single layer operator and double layer
  [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeKernelsOperatorsAxis(Xsing,Ysing,nnn,x,y,PARAM);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %compute influence matrix and known term
  [A,b,nnx,nny] = buildInfluenceMatrixAndKnownTerm(GXX,GXY,GYX,GYY,A11,A12,A21,A22,Xsing,nnn,x,y,PARAM);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
  %deflation, it is needed for inviscid bubbles and closed walls (like particles)
  [A,b] = deflationStokes(Ysing,nnn,x,y,A,b,PARAM);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %add force free condition
  [A,b] = forceFreeStokesAxis(Xsing,Ysing,nnn,x,y,A,b,PARAM);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %add volume constaint like Alinovi and Bottaro
  %[A,b] = volumeCorrectionAlinovi(Xsing,Ysing,nnn,x,y,A,b,PARAM);
  
  %solve linear system
  xxx = A\b; 
  
end
  