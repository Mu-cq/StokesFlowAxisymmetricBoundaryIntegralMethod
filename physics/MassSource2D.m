%determine BC depending also from an eventual mass injection

function [z,ux,uy,fx,fy] = MassSource2D(Xsing,Ysing,PARAM,r,df_x,df_y)

n = PARAM.n;    m = PARAM.m;    j = PARAM.j;    q = PARAM.q;
fixed_elem = n+2*m+j;

ux = 0;
uy = 0;
fx = 0;
fy = 0;

%%%%%%%%%%%%% ADD MASS CONTRIBUTION WITH POTENTIAL FLOW %%%%%%%%%%%%%%%%%
  if PARAM.continuity==1;
  
      ux = 0;
      uy = 0;
      NXX = 0;
      NXY = 0;
      NYX = 0;
      NYY = 0;
      
      for l = 1:numel(PARAM.mass)
  
        %add velocity field and stresses with potential flow
        %[~,u_temp,v_temp,nxx_temp,nxy_temp,nyx_temp,nyy_temp] = gf_2D_Laplace(Xsing,Ysing,PARAM.x_inj(l),PARAM.y_inj(l),...
         %   PARAM.R,-PARAM.R,visc2,PARAM.p_out,PARAM.mass(l));
        [~,ux_temp,uy_temp,nxx_temp,nxy_temp,nyx_temp,nyy_temp] = gf2D_laplace_fs(Xsing,Ysing,PARAM.x_inj(l),PARAM.y_inj(l),PARAM.mass(l));
          
              ux = ux + ux_temp;
              uy = uy + uy_temp;
              NXX = NXX + [nxx_temp(1:fixed_elem) zeros(1,q)];
              NXY = NXY + [nxy_temp(1:fixed_elem) zeros(1,q)];
              NYX = NYX + [nyx_temp(1:fixed_elem) zeros(1,q)];
              NYY = NYY + [nyy_temp(1:fixed_elem) zeros(1,q)];
                         
          %end
    
      end

    %add stresses
    fx = NXX.*r(1,:)+NXY.*r(2,:);
    fy = NYX.*r(1,:)+NYY.*r(2,:);
    
    %vector containing the MODIFIED boundaries conditions
    z=zeros(2*(n+2*m+j+q),1);
    
    z(1:2:2*n-1) = -PARAM.p_out-fx(1:n);   %OUTLET x direction
    z(2:2:2*n) = -uy(1:n);  %OUTLET y direction
    z(2*n+1:2:2*(n+m)-1) = -ux(n+1:n+m);   %UPPER WALL  x direction
    z(2*n+2:2:2*(n+m)) = -uy(n+1:n+m);   %UPPER WALL y direction
    z(2*(n+m)+1:2:2*(n+m+j)-1) = -(PARAM.p_out+PARAM.press_grad)-fx(n+m+1:n+m+j);  %INLET  x direction
    z(2*(n+m)+2:2:2*(n+m+j)) = -uy(n+m+1:n+m+j);  %INLET y direction
    z(2*(n+m+j)+1:2:2*(n+2*m+j)-1) = -ux(n+m+j+1:n+m+m+j);    %LOWER WALL x direction
    z(2*(n+m+j)+2:2:end) = -uy(n+m+j+1:end);    %LOWER WALL y direction
    
    z(2*(n+2*m+j)+1:2:end-1) = df_x - fx(n+2*m+j+1:end);    %INTERFACE X
    z(2*(n+2*m+j)+2:2:end) = df_y - fy(n+2*m+j+1:end);    %INTERFACE Y
    
  else
      %vector containing the PHYSICAL boundaries conditions
      z=zeros(2*(n+2*m+j+q),1);

      %stresses boundary conditions
      z(1:2:2*n-1) = -PARAM.p_out*ones(n,1);
      z(2*(n+m)+1:2:2*(n+m+j)-1) = -PARAM.p_out*ones(n,1)-PARAM.press_grad;
      
      z(2*(n+m+m+j)+1:2:end-1) = df_x';
      z(2*(n+m+m+j)+2:2:end) = df_y';

  end
  
end