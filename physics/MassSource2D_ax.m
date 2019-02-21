%determine BC depending also from an eventual mass injection

function [z,ux,uy,fx,fy] = MassSource2D_ax(Xsing,Ysing,PARAM,r,df_x,df_y)

n = PARAM.n;    m = PARAM.m;    j = PARAM.j;    q = PARAM.q+1;
fixed_elem = n+m+j;

if PARAM.spline == 2
    m = m+1;
    fixed_elem = fixed_elem + 1;
end

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
          
%           hold on
%           plot(PARAM.x_inj(l),PARAM.y_inj(l),'om')
%           hold off
  
        %add velocity field and stresses with potential flow
        [ux_temp,uy_temp,nxx_temp,nxy_temp,nyx_temp,nyy_temp] = gf2D_ax_laplace_fs(Xsing,Ysing,PARAM.x_inj(l),PARAM.mass(l));
          
              ux = ux + ux_temp;
              uy = uy + uy_temp;
              NXX = NXX + PARAM.visc2*[nxx_temp(1:fixed_elem) zeros(1,q)];
              NXY = NXY + PARAM.visc2*[nxy_temp(1:fixed_elem) zeros(1,q)];
              NYX = NYX + PARAM.visc2*[nyx_temp(1:fixed_elem) zeros(1,q)];
              NYY = NYY + PARAM.visc2*[nyy_temp(1:fixed_elem) zeros(1,q)];
                         
          %end
    
      end

    %add stresses
    if PARAM.spline == 2
        fx = NXX.*r(1,:)+NXY.*r(2,:);
        fy = NYX.*r(1,:)+NYY.*r(2,:);
    else
        fx = NXX.*r(1,:)+NXY.*r(2,:);
        fy = NYX.*r(1,:)+NYY.*r(2,:);
    end
    
    %vector containing the MODIFIED boundaries conditions
    z=zeros(2*(n+m+j+q),1);
    
    %if inlet and outlet exist
    if PARAM.stress~=1000
        z(1:2:2*n-1) = -PARAM.p_out-fx(1:n);   %OUTLET x direction
        z(2*(n+m)+1:2:2*(n+m+j)-1) = -(PARAM.p_out+PARAM.press_grad)-fx(n+m+1:n+m+j);  %INLET  x direction
    end
    
    z(2:2:2*n) = -uy(1:n);  %OUTLET y direction
    z(2*n+1:2:2*(n+m)-1) = -ux(n+1:n+m);   %UPPER WALL  x direction
    z(2*n+2:2:2*(n+m)) = -uy(n+1:n+m);   %UPPER WALL y direction
    z(2*(n+m)+2:2:2*(n+m+j)) = -uy(n+m+1:n+m+j);  %INLET y direction
    
    z(2*(n+m+j)+1:2:end-1) = df_x - fx(n+m+j+1:end);    %INTERFACE X
    z(2*(n+m+j)+2:2:end) = df_y - fy(n+m+j+1:end);    %INTERFACE Y
    
  else
      
      if PARAM.stress==1
          %vector containing the PHYSICAL boundaries conditions
          z=zeros(2*(n+m+j+q),1);

          %stresses boundary conditions
          z(1:2:2*n-1) = -PARAM.p_out*ones(n,1);
          z(2*(n+m)+1:2:2*(n+m+j)-1) = -PARAM.p_out*ones(n,1)-PARAM.press_grad;

          z(2*(n+m+j)+1:2:end-1) = df_x';
          z(2*(n+m+j)+2:2:end) = df_y';
      elseif PARAM.stress==0
          %vector containing the PHYSICAL boundaries conditions
          z=zeros(2*(n+m+j+q),1);

          %stresses boundary conditions
          z(1:2:2*n-1) = -PARAM.p_out*ones(n,1);
          z((n+m)*2+1:2:2*(n+m+j)-1) = 2*PARAM.Q/pi/(PARAM.R-PARAM.L/2*sin(PARAM.theta))^2*...
                (1-(Ysing(n+m+1:n+m+j)/(PARAM.R-PARAM.L/2*sin(PARAM.theta))).^2);

          z(2*(n+m+j)+1:2:end-1) = df_x';
          z(2*(n+m+j)+2:2:end) = df_y';
          
      else
          %vector containing the PHYSICAL boundaries conditions
          z=zeros(2*(m+q),1);

          z(2*m+1:2:end-1) = df_x';
          z(2*m+2:2:end) = df_y';
      end

  end
  
end