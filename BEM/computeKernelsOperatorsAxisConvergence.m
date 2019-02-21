%deflation on block where A is the influnce matrix and b is known term
%vector

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeKernelsOperatorsAxisConvergence(Xsing,Ysing,nnn,x,y,PARAM)

  %preallocation
  GXX = cell(numel(PARAM.n),1);
  GXY = cell(numel(PARAM.n),1);
  GYX = cell(numel(PARAM.n),1);
  GYY = cell(numel(PARAM.n),1);
  A11 = cell(numel(PARAM.n),1);
  A12 = cell(numel(PARAM.n),1);
  A21 = cell(numel(PARAM.n),1);
  A22 = cell(numel(PARAM.n),1);
  
  %compute kernels
  if PARAM.cfunction == 0   %compute with MATLAB
      
        for i = 1:numel(PARAM.n)
            
            if i==1
                startMatrix = 1;
            else
                startMatrix = 1+sum(nnn(1:i-1));
            end
            
            if PARAM.orderVariable(i)==0 && PARAM.orderGeometry(i)==0   % integration on constant, straight element
        
                %compute Kernels
                [GXX{i},GXY{i},GYX{i},GYY{i},A11{i},A12{i},A21{i},A22{i}] = computeKernelStokesAxisConst(x{i},y{i},Xsing,Ysing,startMatrix,sum(nnn(1:i)),PARAM.STstokes,PARAM.kernelFreeSpace,PARAM.posWall);
        
            elseif PARAM.orderVariable(i)==1 && PARAM.orderGeometry(i)==0     % integration on linear, straight element

                %compute Kernels
                [GXX{i},GXY{i},GYX{i},GYY{i},A11{i},A12{i},A21{i},A22{i}] = computeKernelStokesAxisLinear(x{i},y{i},Xsing,Ysing,startMatrix,sum(nnn(1:i)),PARAM.STstokes,PARAM.kernelFreeSpace,PARAM.posWall);
                
            elseif PARAM.orderVariable(i)==0 && PARAM.orderGeometry(i)==1     % integration on constant, curved element
                
                %compute spline coeff
                if PARAM.SPlinesType(i)==1
                    [ax,bx,cx,dx,ay,by,cy,dy] = spline_natural(x{i},y{i});
                elseif PARAM.SPlinesType(i)==2
                    [ax,bx,cx,dx,ay,by,cy,dy] = spline_symmetric(x{i},y{i});
                end
                
                error('probably bug in singular treatment')

                %compute Kernels
                [GXX{i},GXY{i},GYX{i},GYY{i},A11{i},A12{i},A21{i},A22{i}] = computeKernelStokesAxisConstSpline(Xsing,Ysing,startMatrix,sum(nnn(1:i)),PARAM.STstokes,ax,ay,bx,by,cx,cy,dx,dy,PARAM.kernelFreeSpace,PARAM.posWall);
                
            elseif PARAM.orderVariable(i)==1 && PARAM.orderGeometry(i)==1     % integration on linear, curved element

                %compute spline coeff
                if PARAM.SPlinesType(i)==1
                    [ax,bx,cx,dx,ay,by,cy,dy] = spline_natural(x{i},y{i});
                elseif PARAM.SPlinesType(i)==2
                    [ax,bx,cx,dx,ay,by,cy,dy] = spline_symmetric(x{i},y{i});
                end
                
                %compute Kernels
                [GXX{i},GXY{i},GYX{i},GYY{i},A11{i},A12{i},A21{i},A22{i}] = computeKernelStokesAxisLinearSplineConvergence(Xsing,Ysing,startMatrix,sum(nnn(1:i)),PARAM.STstokes,ax,ay,bx,by,cx,cy,dx,dy,PARAM.kernelFreeSpace,PARAM.posWall);
                
            end
        
        end


  elseif PARAM.cfunction == 1
          
          error('Not implemented')
      
  end