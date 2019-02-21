%compute normal vector speectrally

function [nx,ny] = normalVectorSpectral(x,y,PARAM)

%derivatives
D1 = PARAM.D1;

%compute geomtrical derivaties
xp = D1*x;    yp = D1*y;
  
%compute normal vector
h = (xp.^2+yp.^2).^(0.5);
nx = yp./h;
ny = -xp./h;
    
end