%visualization 2D Stokes, integration on consatnt elements

function [GXX,GXY,GYY,A11,A12,A22] = GT_2DStokes_LinearSplinesElement(x,y,x0,y0)

      field = numel(x0);
      elem = numel(x)-1;

      %interface normal and curvature
      [ax, bx, cx, dx, ay, by, cy, dy] = my_spline_periodic ([x(1:end-1)' x(1)], [y(1:end-1)' y(1)]);
      
      %integration on interface, linear curved element (SPlines)
      [GXX,GXY,GYY,TXXX,TXXY,TXYY,TYYY] =...
            Stokes2DSPlinesLinear(ax,bx,cx,dx,ay,by,cy,dy,x0',y0');
      
      %reshape
      GXX = reshape(GXX,field,elem);
      GXY = reshape(GXY,field,elem);
      GYY = reshape(GYY,field,elem);
      TXXX = reshape(TXXX,field,elem);
      TXXY = reshape(TXXY,field,elem);
      TXYY = reshape(TXYY,field,elem);
      TYYY = reshape(TYYY,field,elem);
      
      A11 = TXXX+TXXY;
      A12 = TXXY+TXYY;
      A22 = TXYY+TYYY;
      
end