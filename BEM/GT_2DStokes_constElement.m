%visualization 2D Stokes, integration on consatnt elements

function [GXX,GXY,GYY,A11,A12,A22] = GT_2DStokes_constElement(x,y,x0,y0)

      fixed_elem = numel(x)-1;
      field = numel(x0);

      no = sqrt(diff(x).^2+diff(y).^2);
      r1 = -diff(y)./no;
      r2 = diff(x)./no;

      R1 = repmat(r1',field,1);
      R2 = repmat(r2',field,1);
      
      %integration field points on wall
      [GXX,GXY,GYY,TXXX,TXXY,TXYY,TYYY] = Stokes2DGaussIntConst(x',y',x0',y0');
      
      %reshape
      GXX = reshape(GXX,field,fixed_elem);
      GXY = reshape(GXY,field,fixed_elem);
      GYY = reshape(GYY,field,fixed_elem);
      TXXX = reshape(TXXX,field,fixed_elem);
      TXXY = reshape(TXXY,field,fixed_elem);
      TXYY = reshape(TXYY,field,fixed_elem);
      TYYY = reshape(TYYY,field,fixed_elem);
      
      A11 = TXXX.*R1 + TXXY.*R2;
      A12 = TXXY.*R1 + TXYY.*R2;
      A22 = TXYY.*R1 + TYYY.*R2;
      
end