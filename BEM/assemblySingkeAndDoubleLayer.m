%assembly single and double layer operator

function [SL,DL,SLpressure,DLpressure,fBC,uBC] = assemblySingkeAndDoubleLayer(solution,GXX,GXY,GYX,GYY,A11,A12,A21,A22,PX,PY,PI1,PI2,Xsing,nnn,x,y,PARAM)

%initialize
SL = zeros(2*numel(Xsing),2*sum(nnn));
DL = zeros(2*numel(Xsing),2*sum(nnn));
SLpressure = zeros(numel(Xsing),2*sum(nnn));
DLpressure = zeros(numel(Xsing),2*sum(nnn));
uBC = zeros(2*sum(nnn),1);
fBC = zeros(2*sum(nnn),1);

%loop in panels for integration
for i = 1:numel(PARAM.n)
    
      if PARAM.typeBC(i)==2 || PARAM.typeBC(i)==7
          preDL = 1-PARAM.visc(i);
      else
          preDL = 1;
      end
      
      %range
      [startMatrix,endMatrix] = getSingRange(i,nnn);
      
      if (PARAM.typeBCstokes(i)==0||PARAM.typeBCstokes(i)==1||PARAM.typeBCstokes(i)==2||PARAM.typeBCstokes(i)==3||PARAM.typeBCstokes(i)==4||PARAM.typeBCstokes(i)==6||PARAM.typeBCstokes(i)==7)
      
          %single layer
          SL(1:2:end-1,2*startMatrix-1:2:2*endMatrix-1) = GXX{i};
          SL(1:2:end-1,2*startMatrix:2:2*endMatrix) = GXY{i};
          SL(2:2:end,2*startMatrix-1:2:2*endMatrix-1) = GYX{i};
          SL(2:2:end,2*startMatrix:2:2*endMatrix) = GYY{i};

          %double layer
          DL(1:2:end-1,2*startMatrix-1:2:2*endMatrix-1) = A11{i}*preDL;
          DL(1:2:end-1,2*startMatrix:2:2*endMatrix) = A12{i}*preDL;
          DL(2:2:end,2*startMatrix-1:2:2*endMatrix-1) = A21{i}*preDL;
          DL(2:2:end,2*startMatrix:2:2*endMatrix) = A22{i}*preDL;

          %single layer pressure
          SLpressure(:,2*startMatrix-1:2:2*endMatrix-1) = PX{i};
          SLpressure(:,2*startMatrix:2:2*endMatrix) = PY{i};

          %double layer pressure
          DLpressure(:,2*startMatrix-1:2:2*endMatrix-1) = PI1{i}*preDL;
          DLpressure(:,2*startMatrix:2:2*endMatrix) = PI2{i}*preDL;
      
      elseif PARAM.typeBCstokes(i)==5
          
          [nx,ny] = computeNormalVector(x{i},y{i},PARAM.orderVariable(i),PARAM.orderGeometry(i),PARAM.SPlinesType(i));
          cosT = repmat(nx,numel(Xsing),1);
          sinT = repmat(ny,numel(Xsing),1);
          
          %modify single layer for normal and tangent variables
          SL11 = GXX{i}.*cosT + GXY{i}.*sinT;
          SL12 = -GXX{i}.*sinT + GXY{i}.*cosT;
          SL21 = GYX{i}.*cosT + GYY{i}.*sinT;
          SL22 = -GYX{i}.*sinT + GYY{i}.*cosT;

          %double layer for normal and tangent
          DL11 = A11{i}.*cosT + A12{i}.*sinT;
          DL12 = -A11{i}.*sinT + A12{i}.*cosT;
          DL21 = A21{i}.*cosT + A22{i}.*sinT;
          DL22 = -A21{i}.*sinT + A22{i}.*cosT;
          
          %modify PRESSURE single layer for normal and tangent variables
          SL1p = PX{i}.*cosT + PY{i}.*sinT;
          SL2p = -PX{i}.*sinT + PY{i}.*cosT;
          
          %modify PRESSURE double layer for normal and tangent variables
          DL1p = PI1{i}.*cosT + PI2{i}.*sinT;
          DL2p = -PI1{i}.*sinT + PI2{i}.*cosT;
          
          %single layer
          SL(1:2:end-1,2*startMatrix-1:2:2*endMatrix-1) = SL11;
          SL(1:2:end-1,2*startMatrix:2:2*endMatrix) = SL12;
          SL(2:2:end,2*startMatrix-1:2:2*endMatrix-1) = SL21;
          SL(2:2:end,2*startMatrix:2:2*endMatrix) = SL22;

          %double layer
          DL(1:2:end-1,2*startMatrix-1:2:2*endMatrix-1) = DL11*preDL;
          DL(1:2:end-1,2*startMatrix:2:2*endMatrix) = DL12*preDL;
          DL(2:2:end,2*startMatrix-1:2:2*endMatrix-1) = DL21*preDL;
          DL(2:2:end,2*startMatrix:2:2*endMatrix) = DL22*preDL;
          
          %single layer pressure
          SLpressure(:,2*startMatrix-1:2:2*endMatrix-1) = SL1p;
          SLpressure(:,2*startMatrix:2:2*endMatrix) = SL2p;

          %double layer pressure
          DLpressure(:,2*startMatrix-1:2:2*endMatrix-1) = DL1p*preDL;
          DLpressure(:,2*startMatrix:2:2*endMatrix) = DL2p*preDL;
          
      else
          
          error('Not implemented')
          
      end
      
      if PARAM.typeBC(i)==0 || PARAM.typeBC(i)==1 || PARAM.typeBC(i)==3 || PARAM.typeBC(i)==6   % prescribe velocity, stresses are from the result
          
          %compute normal vector
          [nx,ny] = computeNormalVector(x{i},y{i},PARAM.orderVariable(i),PARAM.orderGeometry(i));
          if PARAM.typeBC(i)==0
                velx = PARAM.velBCaxial{i};
                vely = PARAM.velBCradial{i};
          elseif PARAM.typeBC(i)==1
                velx = PARAM.velBC{i}.*nx';
                vely = PARAM.velBC{i}.*ny';
          elseif PARAM.typeBC(i)==3
                tx = -ny;   ty = nx;
                velx = PARAM.velBC{i}.*tx';
                vely = PARAM.velBC{i}.*ty';
          elseif PARAM.typeBC(i)==6
                velx = 0;
                vely = 0;
          end

          %known BC
          uBC(2*startMatrix-1:2:2*endMatrix-1) = velx;
          uBC(2*startMatrix:2:2*endMatrix) = vely;
          
          %BC from the result
          fBC(2*startMatrix-1:2:2*endMatrix-1) = solution(2*startMatrix-1:2:2*endMatrix-1);
          fBC(2*startMatrix:2:2*endMatrix) = solution(2*startMatrix:2:2*endMatrix);
          
      elseif PARAM.typeBC(i)==2 || PARAM.typeBC(i)==7   % prescribe stresses, velocity are from the solution
          
          %compute curvature
          [nx,ny] = computeNormalVector(x{i},y{i},PARAM.orderVariable(i),PARAM.orderGeometry(i),PARAM.SPlinesType(i));
          [k1,k2] = computeCurvatureSplines(x{i},y{i},PARAM.orderVariable(i));
          k = k1+k2;
          
          %known BC
          fBC(2*startMatrix-1:2:2*endMatrix-1) = PARAM.stressBC{i}*k.*nx;
          fBC(2*startMatrix:2:2*endMatrix) = PARAM.stressBC{i}*k.*ny;
          
          %add gravity
          if PARAM.typeBC(i)==7
          
              buoyancy = -x{i}*PARAM.Bond;
              fBC(2*startMatrix-1:2:2*endMatrix-1) = fBC(2*startMatrix-1:2:2*endMatrix-1) + buoyancy'.*nx';
              fBC(2*startMatrix:2:2*endMatrix) = fBC(2*startMatrix:2:2*endMatrix) + buoyancy'.*ny';
              
          end
          
          %BC from the result
          uBC(2*startMatrix-1:2:2*endMatrix-1) = solution(2*startMatrix-1:2:2*endMatrix-1);
          uBC(2*startMatrix:2:2*endMatrix) = solution(2*startMatrix:2:2*endMatrix);
          
      elseif PARAM.typeBC(i)==5    % prescribe normal velocity and tangent stress
          
          %known BC
          fBC(2*startMatrix-1:2:2*endMatrix-1) = solution(2*startMatrix-1:2:2*endMatrix-1);
          fBC(2*startMatrix:2:2*endMatrix) = PARAM.stressBC{i};
          
          %BC from the result
          uBC(2*startMatrix-1:2:2*endMatrix-1) = PARAM.velBC{i};
          uBC(2*startMatrix:2:2*endMatrix) = solution(2*startMatrix:2:2*endMatrix);
          
%       elseif PARAM.typeBC(i)==7    % prescribe normal stress and tangent velocity
%           
%           %known BC
%           fBC(2*startMatrix-1:2:2*endMatrix-1) = solution(2*startMatrix-1:2:2*endMatrix-1);
%           fBC(2*startMatrix:2:2*endMatrix) = PARAM.stressBC{i};
%           
%           %BC from the result
%           uBC(2*startMatrix-1:2:2*endMatrix-1) = PARAM.velBC{i};
%           uBC(2*startMatrix:2:2*endMatrix) = solution(2*startMatrix:2:2*endMatrix);
%           
%       elseif PARAM.typeBC(i)==8    % prescribe normal and tangent stress
%           
%           %known BC
%           fBC(2*startMatrix-1:2:2*endMatrix-1) = solution(2*startMatrix-1:2:2*endMatrix-1);
%           fBC(2*startMatrix:2:2*endMatrix) = PARAM.stressBC{i};
%           
%           %BC from the result
%           uBC(2*startMatrix-1:2:2*endMatrix-1) = PARAM.velBC{i};
%           uBC(2*startMatrix:2:2*endMatrix) = solution(2*startMatrix:2:2*endMatrix);
              
      else
          
          error('Not implemented')
          
      end
      
end