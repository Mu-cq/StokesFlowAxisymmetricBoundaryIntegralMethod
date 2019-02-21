function [nx,ny] = normal_splines2(bx,cx,dx,by,cy,dy)
%compute normals with splines

    %compute the versor normal tp the node (N) and to the element (r)
  %[N,r] =normal_versor_clean(a,b,q-1);
  N = [by./sqrt(bx.*bx+by.*by) (by(end)+2*cy(end)+3*dy(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)));...
      -bx./sqrt(bx.*bx+by.*by) (-bx(end)-2*cx(end)-3*dx(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)))];
  
  nx = N(1,:);  ny = N(2,:);

end

