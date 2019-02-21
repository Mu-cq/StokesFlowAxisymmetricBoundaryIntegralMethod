%give the velocities on a droplet interface when the shape is known

function [sol,nx,ny,Vdrop] = bemChebXY(x,y,PARAM)

  error('Use bemSpectralXY')

  %number of elements
  n = PARAM.n;
  visc = PARAM.visc;
  
  %differentiation matrix
  D1 = PARAM.D1;
  D2 = PARAM.D2;
  
  %integration weight
  WG = PARAM.manyWG;

  %gri points for green function
  X0 = repmat(x,1,n+1);
  Y0 = repmat(y,1,n+1);
  X = X0';
  Y = Y0';
  
  %clean the diagonal for singular treatamet
  Y = Y + eye(n+1);% Y([1 end]) = Y([1 end]) - [1 1];

  %compute green's functions
  [SXX,SXY,SYX,SYY] = axisKernelSingleLayer(X,Y,X0,Y0);
  
  %clean the diagonal for singular treatamet
  SXX = SXX - diag(diag(SXX));
  SXY = SXY - diag(diag(SXY));
  SYX = SYX - diag(diag(SYX));
  SYY = SYY - diag(diag(SYY));
  
  %compute geomtrical derivaties
  xp = D1*x;    yp = D1*y;
  xpp = D2*x;    ypp = D2*y;
  
  %compute normal vector
  h = (xp.^2+yp.^2).^(0.5);
  nx = yp./h;
  ny = -xp./h;
  
  %compute curvature in meridional plane
  K1 = (xp.*ypp-yp.*xpp)./(xp.^2+yp.^2).^(1.5);
  
  %compute curvature in aximuthal direction
  K2 = ny./y;
  K2([1 end]) = K1([1 end]);
  
  %total curvature
  K = K1+K2;
  
  if PARAM.BC==2
      
      cap_forces = K/PARAM.Ca;
      buoyancy = x*3*(1+1.5*visc)/(1+visc);
      
      K = cap_forces + buoyancy;
      
  end
  
  %modify curvature for singualr treatment (CURVATURE, NOT STRESSES!!!)
  prepareK = repmat(K',2,1);
  prepareK2 = repmat(prepareK(:),1,2);
  KKK = repmat(prepareK2',n+1,1)-repmat(prepareK2,1,n+1);
  
  %in this case gamma=1 by definition, as Stone et al. 1989
  KKK = KKK*1;
  
  %prepare for integration
  NX = repmat(nx',2*n+2,1);
  NY = repmat(ny',2*n+2,1);
  prepareH = repmat(h',2,1);
  H = repmat(prepareH(:)',2*n+2,1);
  
  %RHS integration
  U = zeros(2*(n+1));
  zzz = zeros(2*(n+1));
  zzz(:,1:2:end-1) = NX;
  zzz(:,2:2:end) = NY;
  U(1:2:end-1,1:2:end-1) = SXX;
  U(1:2:end-1,2:2:end) = SXY;
  U(2:2:end,1:2:end-1) = SYX;
  U(2:2:end,2:2:end) = SYY;
  b = -sum(U.*zzz.*KKK.*WG.*H,2);   %integration using formual in reasearch notes
  
  if PARAM.BC==1
      %add underlying flow
      [ux,uy] = extens_flow(x',y',PARAM.Ca);
      b(1:2:end-1) = b(1:2:end-1) + 8*pi*ux;
      b(2:2:end) = b(2:2:end) + 8*pi*uy;
  end
  %LHS matrix
  %A = zeros(2*(n+1));
  %A = A + 4*pi*(1+visc)*eye(2*n+2);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %find the solution in terms of stresses and velocities with matrix
  %inversion Ax = b
  
  %solve linear system
  %x = A\b;
  sol = b./(4*pi*(1+visc)*ones(2*n+2,1));
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %compute droplet velocity
  ux = sol(1:2:end-1);    uy = sol(2:2:end);
  uNormal = ux.*nx + uy.*ny;
  Vdrop = DropVelocityAxisCurvilinear(x,y,uNormal,PARAM);
  
  if PARAM.dropFrame==1
      sol(1:2:end-1) = sol(1:2:end-1)-Vdrop;
  end
  
  %T = toc;
  
end
  