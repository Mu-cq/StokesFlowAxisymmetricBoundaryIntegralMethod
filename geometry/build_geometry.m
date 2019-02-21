%build geometry

function [a,b,MIN] = build_geometry2(ampli,D,q,geometry,lealpoz,visc,Ca)

    if (geometry==1)

          %%%%%%%%%%%%%%%%% BUILD DOMAIN %%%%%%%%%%%%%%%%%%
          %build boundary with perturbation
          if (lealpoz==1)
              [a,b,~,MIN] = draw_circle_lean(0,0,q,D);
          elseif (lealpoz==0)
              [a,b,~,MIN] = draw_circle_clean(0,0,q,ampli);
          end
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      elseif (geometry==2)

          %%%%%%%%%%%%%%%%% IMPORT DOMAIN CONFIGURATION %%%%%%%%%%%
          cd domain
          load savea.mat
          load saveb.mat
          load elem
          load ds.mat
          cd ..

          MIN = ds;

          if (lealpoz==1)
              D = 1000;
          elseif (lealpoz==0)
              ampli = 1000;
          end
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      elseif (geometry==3)

          %%%%%%%%%%%%%%%%% IMPORT OPTIMAL SHAPE %%%%%%%%%%%%%%%%%%
          cd optimal_shape
          folder = strcat('lambda=',num2str(visc),'_Ca=',num2str(Ca));
          cd(folder)

          filename = strcat('crd_rzlam',num2str(visc),'_ca',num2str(Ca),'_mA1000_delta0.02_time0.dat');

          load(filename)
          x = crd_rzlam0_5_ca12_mA1000_delta0_02_time0(:,1);
          y = crd_rzlam0_5_ca12_mA1000_delta0_02_time0(:,2);

          cd ..
          cd ..

          %[a, b] = remesh7(y,x,q);
          [a, b] = adapt(y,x,q);

          %compute MIN
          %ds = sqrt((a(1:end-1)-a(2:end)).^2+(b(1:end-1)-b(2:end)).^2);

          [xxx, yyy] = remesh7(y,x,q);
          ds = sqrt((xxx(1:end-1)-xxx(2:end)).^2+(yyy(1:end-1)-yyy(2:end)).^2);

          MIN=min(ds);
          D = 1000;
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
  
end