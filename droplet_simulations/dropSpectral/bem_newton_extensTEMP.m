%give the solution of the poiseiulle flow in a pipe in terms of stress and velocities
%having the following input: radius of the pipe R, length of the pipe L,
%number of element at the outlet n, number of element at the wall m, number of
%element at the inlet j, inlet velocity
%vel_in, vertical stress at the outlet "stress" and the viscosity of the
%fluid



function [y,N,ax,bx,cx,dx,ay,by,cy,dy,df_x,df_y,K1,K2]=bem_newton_extens(a,b,q,lambda,capillary)

  %dimensional viscosity
  visc2 = 1;
  visc1 = lambda;

  %from elements to nodes
  q = q+1;

  %tic
  
  %compute the spline coeff
  [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (a, b);
  
  N = [by./sqrt(bx.*bx+by.*by) by(end)+2*cy(end)+3*dy(end)/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)));...
      -bx./sqrt(bx.*bx+by.*by) -bx(end)-2*cx(end)-3*dx(end)/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)))];
  
  %%%%%%%%%%%%%%%%%%% COMPUTE CURVATURE WITH SPLINES %%%%%%%%%%%%%%%%%%%%%%

  K1 = curv_spline2(bx,by,cx,cy,dx,dy);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %those because the points are on the axis
  N(:,1) = [1; 0];
  N(:,end) = [-1; 0];
  
  %second component of the curvature
  K2 = N(2,:)./b;
  K2(1) = K1(1);
  K2(end) = K1(end);
  
  %sum the two component of the curvature
  K = K1+K2;
  
  %calculate deltaf
  df = K;
  [df_x,df_y] = stress_diff(df,N,q);
  
  %compute Gree's function
  %if PARAM.cfunction == 0

      %[GXX,GXY,GYX,GYY,A11,A12,A21,A22,T1,T2,D1,D2] = computeGT_spline_vect(a,b,ax,ay,bx,by,cx,cy,dx,dy);

  %elseif PARAM.cfunction == 1
           

      [GXX,GXY,GYX,GYY,TXXX,TXXY,TXYX,TXYY,TYXX,TYXY,TYYX,TYYY,T1,T2,D1,D2] = ...
              Stokes2DAxisSPlinesLinear(ax,bx,cx,dx,ay,by,cy,dy,a,b);
              
      %reshape because c function gives a vector
      GXX = reshape(GXX,q,q);
      GXY = reshape(GXY,q,q);
      GYX = reshape(GYX,q,q);
      GYY = reshape(GYY,q,q);
      %plus because normal is pointing out
      TXXX = reshape(TXXX,q,q);
      TXXY = reshape(TXXY,q,q);
      TXYX = reshape(TXYX,q,q);
      TXYY = reshape(TXYY,q,q);
      TYXX = reshape(TYXX,q,q);
      TYXY = reshape(TYXY,q,q);
      TYYX = reshape(TYYX,q,q);
      TYYY = reshape(TYYY,q,q);
      T1 = reshape(T1,q,q-1);
      T2 = reshape(T2,q,q-1);
      D1 = reshape(D1,q,q-1);
      D2 = reshape(D2,q,q-1);

      A11 = TXXX + TXXY;
      A12 = TXYX + TXYY;
      A21 = TYXX + TYXY;
      A22 = TYYX + TYYY;
      
  %end
  
  [u_x,u_y] = extens_flow(a,b,capillary);
  %force symmetry conditions
  u_y(1) = 0;
  u_y(end) = 0;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %compute the vector v containing the known terms:
  
  %vector containing the boundaries conditions
  z=zeros(2*q,1);

  z(1:2:end-1) = df_x;
  z(2:2:end) = df_y;
  
  %preallocation
  U=zeros(2*q);
  
  %right hand side
  U(1:2:end-1,1:2:end-1) = -GXX;
  U(1:2:end-1,2:2:end) = -GXY;
  U(2:2:end,1:2:end-1) = -GYX;
  U(2:2:end,2:2:end) = -GYY;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %vector containing the known terms
  v = -U*z;
  v(1:2:end-1) = v(1:2:end-1) - 8*pi*u_x;
  v(2:2:end) = v(2:2:end) - 8*pi*u_y;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %compute the A matrix containing the coefficients of the unknowns of the
  %problem  
  
  A = zeros(2*q);
  
  %compose A12 and A22  
  A(1:2:end-1,1:2:end-1) = (visc2-visc1)*A11 + diag((-sum(T1,2)-4*pi)*(visc2-visc1) - 4*pi*(visc2+visc1));
  A(1:2:end-1,2:2:end) = (visc2-visc1)*A12 + diag(-sum(D1,2)*(visc2-visc1));
  A(2:2:end,1:2:end-1) = (visc2-visc1)*A21 + diag((-sum(T2,2))*(visc2-visc1));
  A(2:2:end,2:2:end) = (visc2-visc1)*A22 + diag((-sum(D2,2)-4*pi)*(visc2-visc1) - 4*pi*(visc2+visc1));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %find the solution in terms of stresses and velocities with matrix
  %inversion
  
  %v(end) = 0;
  y = A\v;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %T = toc;
  
end
  