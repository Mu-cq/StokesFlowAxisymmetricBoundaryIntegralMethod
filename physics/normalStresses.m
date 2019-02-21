%compute normal stresses

function [dfX,dfY] = normalStresses(x,y,gamma)

[~, bx, cx, dx, ~, by, cy, dy] = spline_symmetric(x,y);
[nx,ny] = normal_splines2(bx,cx,dx,by,cy,dy);
K1 = curv_spline2(bx,by,cx,cy,dx,dy);
K2 = ny./y; K2([1 end]) = K1([1 end]);
dfX = nx.*(K1+K2)*gamma;
dfY = ny.*(K1+K2)*gamma;