%give the solution in terms of stress and velocities
%having the following input: radius of the pipe R, length of the pipe L,
%number of element at the outlet n, number of element at the wall m, number of
%element at the inlet j, inlet velocity
%vel_in, vertical stress at the outlet "stress" and the viscosity of the
%fluid

function [x,Xsing,Ysing,nnn] = BEM_Laplace(x,y,PARAM)

  %change variables name
  PARAM.typeBC = PARAM.typeBClaplace;
  PARAM.orderVariable = PARAM.orderVariableLaplace;
  PARAM.orderGeometry = PARAM.orderGeometryLaplace;
  
  %compute location of the singularity
  [Xsing,Ysing,nnn] = computeSingularityLocation(x,y,PARAM);
  
  %preallocation
  G1 = cell(numel(PARAM.n),1);
  G2 = cell(numel(PARAM.n),1);
  
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
                [G1{i},G2{i}] = computeKernelLaplaceAxis2(x{i},y{i},Xsing,Ysing,startMatrix,sum(nnn(1:i)),PARAM.STlaplace);
        
            elseif PARAM.orderVariable==1 && PARAM.orderGeometry==0     % integration on linear, straight element

                %compute Kernels
                [G1{i},G2{i}] = computeKernelLaplaceAxisLinear(x{i},y{i},Xsing,Ysing,startMatrix,sum(nnn(1:i)),PARAM.STlaplace);

            elseif PARAM.orderVariable==0 && PARAM.orderGeometry==1     % integration on constant, curved element

                error('Not implemented')

            elseif PARAM.orderVariable==1 && PARAM.orderGeometry==1     % integration on linear, curved element

                error('Not implemented')

            end
        
        end


  elseif PARAM.cfunction == 1
          
          %my c function, needs the coordinates of the elements and the field
          %coordiantes (for const elements is the middle point)
%           [G1,Gx,Gy] = Poisson2DConst(a,b,(a(1:end-1)+a(2:end))/2,(b(1:end-1)+b(2:end))/2);
% 
%           G1 = reshape(G1,fixed_elem,fixed_elem);
%           Gx = reshape(Gx,fixed_elem,fixed_elem);
%           Gy = reshape(Gy,fixed_elem,fixed_elem);
      
  end

%   R1 = repmat(r(1,1:end),fixed_elem,1);
%   R2 = repmat(r(2,1:end),fixed_elem,1);
% 
%   G2 = Gx.*R1 + Gy.*R2;
  

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %compute the vector v containing the known terms:
  
  %preallocation
  %U = zeros(fixed_elem);
  
  %build matrices, solve Ax + Uv = 0 where x are the unkwons and v the
  %imposed BC
  A = zeros(numel(Xsing));
  U = zeros(numel(Xsing));
  v = zeros(numel(Xsing),1);
  for i = 1:numel(PARAM.n)
      
      if i==1
          startMatrix = 1;
      else
          startMatrix = 1+sum(nnn(1:i-1));
      end
      
      if PARAM.typeBC(i)==1    %prescribe phi and find gradient
          
          %double layer and single layer potential
          U(:,startMatrix:sum(nnn(1:i))) = G2{i};
          A(:,startMatrix:sum(nnn(1:i))) = -G1{i};
          
          %add term on the diagonal
          U(startMatrix:sum(nnn(1:i)),startMatrix:sum(nnn(1:i))) = U(startMatrix:sum(nnn(1:i)),startMatrix:sum(nnn(1:i))) - diag(0.5*ones(nnn(i),1));
          
          %known BC
          v(startMatrix:sum(nnn(1:i))) = ones(nnn(i),1).*PARAM.concBC{i};
          
      elseif PARAM.typeBC(i)==2    %prescribe gradient and find phi
          
          %single layer and double layer potential
          U(:,startMatrix:sum(nnn(1:i))) = -G1{i};
          A(:,startMatrix:sum(nnn(1:i))) = G2{i};
          
          %add term on the diagonal
          A(startMatrix:sum(nnn(1:i)),startMatrix:sum(nnn(1:i))) = A(startMatrix:sum(nnn(1:i)),startMatrix:sum(nnn(1:i))) - diag(0.5*ones(nnn(i),1));
          
          %known BC
          v(startMatrix:sum(nnn(1:i))) = ones(nnn(i),1).*PARAM.fluxBC{i};
          
      end
      
  end
  
  %solve linear system Ay = b, where b = -Uv
  b = -U*v;
  x = A\b;
  
end
  