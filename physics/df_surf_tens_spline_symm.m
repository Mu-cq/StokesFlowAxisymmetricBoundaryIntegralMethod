%compute jump of nornal stresses at the interface when having surface
%tension and buoyancy using spline for the normals

function [df_x,df_y,K,K1,K2] = df_surf_tens_spline_symm(a,b,gamma)

    %compute the spline coeff
    [~, bx, cx, dx, ~, by, cy, dy] = spline_symmetric (a, b);
    
    %compute the versor normal tp the node (N) and to the element (r)
    %[N,r] =normal_versor_clean(a,b,q-1);
    N = normal_splines(bx,cx,dx,by,cy,dy);
  
    %those because the points are on the axis
    N(:,1) = [1; 0];
    N(:,end) = [-1; 0];
  
    %%%%%%%%%%%%%% COMPUTE MERIDIONAL CURVATURE WITH SPLINES %%%%%%%%%%%%%%%%

    K1 = curv_spline2(bx,by,cx,cy,dx,dy);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %%%%%%%%%%%%%%%%%%%% COMPUTE AZIMUTHAL CURVATURE %%%%%%%%%%%%%%%%%%%%%%%%
  
    %second component of the curvature
    K2 = N(2,:)./b;
  
    %this because on the axis the I cannot use the definition of before
    K2(1) = K1(1);
    K2(end) = K1(end);
  
    %sum the two component of the curvature
    K = K1+K2;
    
    %BC for stresse at the interface
    cap_forces = K*gamma;
  
    df = cap_forces;
    [df_x,df_y] = stress_diff(df,N,numel(a));

end