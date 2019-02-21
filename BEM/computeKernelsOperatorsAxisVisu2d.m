%deflation on block where A is the influnce matrix and b is known term
%vector

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22,PX,PY,PI1,PI2] = computeKernelsOperatorsVisu2d(Xsing,Ysing,x,y,PARAM)

%preallocation
GXX = cell(numel(PARAM.n),1);
GXY = cell(numel(PARAM.n),1);
GYX = cell(numel(PARAM.n),1);
GYY = cell(numel(PARAM.n),1);
A11 = cell(numel(PARAM.n),1);
A12 = cell(numel(PARAM.n),1);
A21 = cell(numel(PARAM.n),1);
A22 = cell(numel(PARAM.n),1);
PX = cell(numel(PARAM.n),1);
PY = cell(numel(PARAM.n),1);
PI1 = cell(numel(PARAM.n),1);
PI2 = cell(numel(PARAM.n),1);
for i = 1:numel(PARAM.n)
            
            if PARAM.orderVariable(i)==0 && PARAM.orderGeometry(i)==0   % integration on constant, straight element
        
                %compute with MATLAB
                [GXX{i},GXY{i},GYX{i},GYY{i},A11{i},A12{i},A21{i},A22{i}] = computeKernelStokesAxisConstNoSing2d(x{i},y{i},Xsing,Ysing);
                [PX{i},PY{i},PI1{i},PI2{i}] = computeKernelPressureStokesConst2d(x{i},y{i},Xsing,Ysing);
        
            elseif PARAM.orderVariable(i)==1 && PARAM.orderGeometry(i)==0     % integration on linear, straight element

                error('not implemented')
                
                [GXX{i},GXY{i},GYX{i},GYY{i},A11{i},A12{i},A21{i},A22{i}] = computeKernelStokesAxisLinearNoSing(x{i},y{i},Xsing,Ysing);
                [PX{i},PY{i},PI1{i},PI2{i}] = computeKernelPressureStokesAxisLinear(x{i},y{i},Xsing,Ysing);

            elseif PARAM.orderVariable(i)==0 && PARAM.orderGeometry(i)==1     % integration on constant, curved element

                error('not implemented')
                
                %compute spline coeff
                if PARAM.SPlinesType(i)==1
                      [ax,bx,cx,dx,ay,by,cy,dy] = spline_natural(x{i},y{i});
                elseif PARAM.SPlinesType(i)==2
                      [ax,bx,cx,dx,ay,by,cy,dy] = spline_symmetric(x{i},y{i});
                end
                
                %compute kernels
                [GXX{i},GXY{i},GYX{i},GYY{i},A11{i},A12{i},A21{i},A22{i}] = computeKernelStokesAxisConstSplineNoSing(Xsing,Ysing,ax,ay,bx,by,cx,cy,dx,dy);
                [PX{i},PY{i},PI1{i},PI2{i}] = computeKernelPressureStokesAxisConstSplineNoSing(Xsing,Ysing,ax,ay,bx,by,cx,cy,dx,dy);

            elseif PARAM.orderVariable(i)==1 && PARAM.orderGeometry(i)==1     % integration on linear, curved element

                error('not implemented')
                
                %compute spline coeff
                if PARAM.SPlinesType(i)==1
                      [ax,bx,cx,dx,ay,by,cy,dy] = spline_natural(x{i},y{i});
                elseif PARAM.SPlinesType(i)==2
                      [ax,bx,cx,dx,ay,by,cy,dy] = spline_symmetric(x{i},y{i});
                end
                
                %compute kernels
                [GXX{i},GXY{i},GYX{i},GYY{i},A11{i},A12{i},A21{i},A22{i}] = computeKernelStokesAxisLinearSplineNoSing(Xsing,Ysing,ax,ay,bx,by,cx,cy,dx,dy);
                [PX{i},PY{i},PI1{i},PI2{i}] = computeKernelPressureStokesAxisLinearSplineNoSing(Xsing,Ysing,ax,ay,bx,by,cx,cy,dx,dy);

            end
        
end