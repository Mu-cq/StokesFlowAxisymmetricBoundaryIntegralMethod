%give the velocities on a droplet interface when the shape is known

function [sol,nx,ny,res,K1,K2,Vdrop,K] = bemSpectralXY(x,y,PARAM)

  %number of elements
  n = PARAM.n;
  visc = PARAM.visc;
  
  %differentiation matrix
  D1 = PARAM.D1;
  D2 = PARAM.D2;
  
  %integration weight
  WG = PARAM.manyWG;

  %grid points for green function
  X0 = repmat(x,1,n+1);
  Y0 = repmat(y,1,n+1);
  X = X0';
  Y = Y0';

  if visc==1    %isoviscosity
      
    %clean the diagonal for singular treatamet
    Y = Y + eye(n+1);
      
    %compute green's functions
    [SXX,SXY,SYX,SYY] = axisKernelSingleLayer(X,Y,X0,Y0);
      
    %clean the diagonal for singular treatamet
    SXX = SXX - diag(diag(SXX));
    SXY = SXY - diag(diag(SXY));
    SYX = SYX - diag(diag(SYX));
    SYY = SYY - diag(diag(SYY));
      
  else
      
    %clean the diagonal for singular treatamet
    Y = Y + eye(n+1);
      
    %compute green's functions
    [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY] = sgf_ax_fs_vect3 (X,Y,X0,Y0);
  
  end
  
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
  
  %in case point goes on the interface
  if PARAM.legendre==2||PARAM.legendre==0
    K2([1 end]) = K1([1 end]);
  end
  
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
  b = -sum(U.*zzz.*KKK.*WG.*H,2);   %integration using formula in reasearch notes
  
  if PARAM.BC==1 || PARAM.BC==3
      %add underlying flow
      if PARAM.BC==1    %linear extensional flow
          [ux,uy] = extens_flow(x',y',PARAM.Ca);
      elseif PARAM.BC==3    %linear extensional flow
          [ux,uy] = nonLinearExtensFlow(x,y,PARAM.Ca,PARAM.CaNL);
      end
      b(1:2:end-1) = b(1:2:end-1) + 8*pi*ux;
      b(2:2:end) = b(2:2:end) + 8*pi*uy;
  end
  
  if visc==1    %isoviscosity
      
      %solution, no need of matrix inversion for isoviscosity
      sol = b./(4*pi*(1+visc)*ones(2*n+2,1));
      
  else
     
     %LHS matrix (double layer for one drop alone)
      A = zeros(2*(n+1));

      %many normal vector
      nnx = repmat(nx',numel(nx),1);
      nny = repmat(ny',numel(ny),1);

      %double layer integration
      DL11 = QXXX.*nnx + QXXY.*nny;
      DL12 = QXYX.*nnx + QXYY.*nny;
      DL21 = QYXX.*nnx + QYXY.*nny;
      DL22 = QYYX.*nnx + QYYY.*nny;

      %from the generalized expansion
      rp = -2./K1;
      
      %store value on the axis
      if PARAM.legendre==0||PARAM.legendre==2
          
          store11 = [DL11(1,1) DL11(end,end)];                        %store value on the axis
          store12 = [DL12(1,1) DL12(end,end)];                        %store value on the axis
          store21 = [DL21(1,1) DL21(end,end)];                        %store value on the axis
          store22 = [DL22(1,1) DL22(end,end)];                        %store value on the axis
          
      end

      %assign values on the diagonal from generalized expansion
      point = 2*ny.^2.*(4*diag(Y0)-ny.*rp)./diag(Y0)./rp;
      DL11(logical(eye(n+1))) = point;

      point = -2*nx.*ny.*(4*diag(Y0)-ny.*rp)./diag(Y0)./rp;
      DL12(logical(eye(n+1))) = point;

      point = -2*nx.*ny.*(4*diag(Y0)-ny.*rp)./diag(Y0)./rp;
      DL21(logical(eye(n+1))) = point;

      point = (16*y.*nx.^2+rp.*(11*ny+ny.^3-3*nx.^2.*ny))./2./diag(Y0)./rp;
      DL22(logical(eye(n+1))) = point;

      %put previous values on the axis
      if PARAM.legendre==0||PARAM.legendre==2
          
         DL11(1,1) = store11(1);   DL11(end,end) = store11(2);
         DL12(1,1) = store12(1);   DL12(end,end) = store12(2);
         DL21(1,1) = store21(1);   DL21(end,end) = store21(2);
         DL22(1,1) = store22(1);   DL22(end,end) = store22(2);
         
      end
      
      %build LHS matrix
      A(1:2:end-1,1:2:end-1) = DL11;
      A(1:2:end-1,2:2:end) = DL12;
      A(2:2:end,1:2:end-1) = DL21;
      A(2:2:end,2:2:end) = DL22;

      %integration weights and metrics term
      A = (1-visc)*A.*WG.*H;

      %add 4*pi*(1+visc)
      A = 4*pi*(1+visc)*eye(2*n+2) - A;
      
      %deflation like Pozrikidis
      if visc<0.1

          %normal vector
          nnn = zzz;

          %operator for integration
          prepareY = repmat(y',2,1);
          yyy = repmat(prepareY(:)',2*n+2,1);
          INTop = 2*pi*yyy.*WG.*H.*nnn'.*nnn;
          A = A - INTop;
          
          %impose volume groth
          b = b - PARAM.Qdeflation*diag(nnn);

      end
      
      %solve linear system
      sol = A\b;
      
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %compute droplet velocity
  ux = sol(1:2:end-1);    uy = sol(2:2:end);
  uNormal = ux.*nx + uy.*ny;
  Vdrop = DropVelocityAxisCurvilinear(x,y,uNormal,PARAM);
  
  %compute residuals
  res = (sol(1:2:end-1)-Vdrop).*nx + sol(2:2:end).*ny;
  res = max(abs(res));
  
  if PARAM.dropFrame==1 || PARAM.dropFrame==2
      sol(1:2:end-1) = sol(1:2:end-1)-Vdrop;
  end
  
end
  