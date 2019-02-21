%deflation on block where A is the influnce matrix and b is known term
%vector

function [A,b,nnx,nny] = buildInfluenceMatrixAndKnownTerm(GXX,GXY,GYX,GYY,A11,A12,A21,A22,Xsing,nnn,x,y,PARAM)

  %initialize
  nnx = cell(numel(PARAM.n),1);
  nny = cell(numel(PARAM.n),1);

  %build matrices, solve Ax + Uv = 0 where x are the unkwons and v the
  %imposed BC
  addSize = sum(PARAM.blockType==1);
  A = zeros(2*numel(Xsing)+addSize);
  U = zeros(2*numel(Xsing)+addSize);
  v = zeros(2*numel(Xsing)+addSize,1);
  for i = 1:numel(PARAM.n)
      
      if i==1
          startMatrix = 1;
      else
          startMatrix = 1+sum(nnn(1:i-1));
      end
      endMatrix = sum(nnn(1:i));
      
      if PARAM.typeBC(i)==0||PARAM.typeBC(i)==1||PARAM.typeBC(i)==3||PARAM.typeBC(i)==4||PARAM.typeBC(i)==6    % prescribe velocity
          
          %double layer and single layer potential
          U(1:2:end-1-addSize,2*startMatrix-1:2:2*endMatrix-1) = A11{i};
          U(1:2:end-1-addSize,2*startMatrix:2:2*endMatrix) = A12{i};
          U(2:2:end-addSize,2*startMatrix-1:2:2*endMatrix-1) = A21{i};
          U(2:2:end-addSize,2*startMatrix:2:2*endMatrix) = A22{i};
          A(1:2:end-1-addSize,2*startMatrix-1:2:2*endMatrix-1) = -GXX{i};
          A(1:2:end-1-addSize,2*startMatrix:2:2*endMatrix) = -GXY{i};
          A(2:2:end-addSize,2*startMatrix-1:2:2*endMatrix-1) = -GYX{i};
          A(2:2:end-addSize,2*startMatrix:2:2*endMatrix) = -GYY{i};
          
          %add term on the diagonal
          if PARAM.panelType(i)==0 || PARAM.panelType(i)==1 % wall
              
                U(2*startMatrix-1:2*endMatrix,2*startMatrix-1:2*endMatrix) = U(2*startMatrix-1:2*endMatrix,2*startMatrix-1:2*endMatrix) - diag(4*pi*ones(2*nnn(i),1));
                
          end
          
          %known BC
          if PARAM.typeBC(i)==0        % axial and radial velocity

                %impose normal velocity
                v(2*startMatrix-1:2:2*endMatrix-1) = PARAM.velBCaxial{i};
                v(2*startMatrix:2:2*endMatrix) = PARAM.velBCradial{i};
                
          elseif PARAM.typeBC(i)==1        % normal velocity
              
                %compute normal vector
                [nx,ny] = computeNormalVector(x{i},y{i},PARAM.orderVariable(i),PARAM.orderGeometry(i),PARAM.SPlinesType(i));

                %impose normal velocity
                v(2*startMatrix-1:2:2*endMatrix-1) = PARAM.velBC{i}.*nx';
                v(2*startMatrix:2:2*endMatrix) = PARAM.velBC{i}.*ny';
                
          elseif PARAM.typeBC(i)==3    % tangent velocity
              
                %compute normal vector
                [nx,ny] = computeNormalVector(x{i},y{i},PARAM.orderVariable(i),PARAM.orderGeometry(i),PARAM.SPlinesType(i));
                tx = -ny;   ty = nx;

                %impose tangent velocity
                v(2*startMatrix-1:2:2*endMatrix-1) = PARAM.velBC{i}.*tx';
                v(2*startMatrix:2:2*endMatrix) = PARAM.velBC{i}.*ty';
                
          elseif PARAM.typeBC(i)==4    % axial velocity

                %impose axial velocity
                v(2*startMatrix-1:2:2*endMatrix-1) = PARAM.velBCaxial{i};
                v(2*startMatrix:2:2*endMatrix) = 0;
                
          end
          
      elseif PARAM.typeBC(i)==2 || PARAM.typeBC(i)==7    %prescribe normal stresses due to surface tension
          
          %error('not implemented')
          visc = PARAM.visc(i);
          
          %compute curvature
          [nx,ny] = computeNormalVector(x{i},y{i},PARAM.orderVariable(i),PARAM.orderGeometry(i),PARAM.SPlinesType(i));
          [k1,k2] = computeCurvatureSplines(x{i},y{i},PARAM.orderVariable(i));
          k = k1+k2;
          
          %double layer and single layer potential
          A(1:2:end-1-addSize,2*startMatrix-1:2:2*endMatrix-1) = (1-visc)*A11{i};
          A(1:2:end-1-addSize,2*startMatrix:2:2*endMatrix) = (1-visc)*A12{i};
          A(2:2:end-addSize,2*startMatrix-1:2:2*endMatrix-1) = (1-visc)*A21{i};
          A(2:2:end-addSize,2*startMatrix:2:2*endMatrix) = (1-visc)*A22{i};
          U(1:2:end-1-addSize,2*startMatrix-1:2:2*endMatrix-1) = -GXX{i};
          U(1:2:end-1-addSize,2*startMatrix:2:2*endMatrix) = -GXY{i};
          U(2:2:end-addSize,2*startMatrix-1:2:2*endMatrix-1) = -GYX{i};
          U(2:2:end-addSize,2*startMatrix:2:2*endMatrix) = -GYY{i};
          
          %add term on the diagonal
          A(2*startMatrix-1:2*endMatrix,2*startMatrix-1:2*endMatrix) = A(2*startMatrix-1:2*endMatrix,2*startMatrix-1:2*endMatrix) - diag(4*pi*ones(2*nnn(i),1))*(1+visc);
          
          %known BC
          v(2*startMatrix-1:2:2*endMatrix-1) = PARAM.stressBC{i}*k.*nx;
          v(2*startMatrix:2:2*endMatrix) = PARAM.stressBC{i}*k.*ny;
          
          %add gravity
          if PARAM.typeBC(i)==7
          
              buoyancy = -x{i}*PARAM.Bond;
              v(2*startMatrix-1:2:2*endMatrix-1) = v(2*startMatrix-1:2:2*endMatrix-1) + buoyancy'.*nx';
              v(2*startMatrix:2:2*endMatrix) = v(2*startMatrix:2:2*endMatrix) + buoyancy'.*ny';
              
          end
          
          %save normal vector
          nnx{i} = nx;
          nny{i} = ny;
          
      elseif PARAM.typeBC(i)==5    %prescribe normal velocity and tangent stress

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
          
          %double layer and single layer potential
          A(1:2:end-1-addSize,2*startMatrix-1:2:2*endMatrix-1) = -SL11;
          A(1:2:end-1-addSize,2*startMatrix:2:2*endMatrix) = DL12;
          A(2:2:end-addSize,2*startMatrix-1:2:2*endMatrix-1) = -SL21;
          A(2:2:end-addSize,2*startMatrix:2:2*endMatrix) = DL22;
          U(1:2:end-1-addSize,2*startMatrix-1:2:2*endMatrix-1) = DL11;
          U(1:2:end-1-addSize,2*startMatrix:2:2*endMatrix) = -SL12;
          U(2:2:end-addSize,2*startMatrix-1:2:2*endMatrix-1) = DL21;
          U(2:2:end-addSize,2*startMatrix:2:2*endMatrix) = -SL22;
          
          %add singular contribution
          A(2*startMatrix-1:2:2*endMatrix-1,2*startMatrix:2:2*endMatrix) = A(2*startMatrix-1:2:2*endMatrix-1,2*startMatrix:2:2*endMatrix) + diag(4*pi*ones(nnn(i),1)).*sinT(startMatrix:endMatrix,:);
          A(2*startMatrix:2:2*endMatrix,2*startMatrix:2:2*endMatrix) = A(2*startMatrix:2:2*endMatrix,2*startMatrix:2:2*endMatrix) - diag(4*pi*ones(nnn(i),1)).*cosT(startMatrix:endMatrix,:);
          U(2*startMatrix-1:2:2*endMatrix-1,2*startMatrix-1:2:2*endMatrix-1) = U(2*startMatrix-1:2:2*endMatrix-1,2*startMatrix-1:2:2*endMatrix-1) - diag(4*pi*ones(nnn(i),1)).*cosT(startMatrix:endMatrix,:);
          U(2*startMatrix:2:2*endMatrix,2*startMatrix-1:2:2*endMatrix-1) = U(2*startMatrix:2:2*endMatrix,2*startMatrix-1:2:2*endMatrix-1) - diag(4*pi*ones(nnn(i),1)).*sinT(startMatrix:endMatrix,:);
          
          %known BC
          v(2*startMatrix-1:2:2*endMatrix-1) = PARAM.velBC{i};
          v(2*startMatrix:2:2*endMatrix) = PARAM.stressBC{i};
          
          %save normal vector
          nnx{i} = nx;
          nny{i} = ny;
          
      else
          
          error('Not implemented')
          
      end
      
  end
  
  %rigid body motion
  under = zeros(2*numel(Xsing)+addSize,1);
  if PARAM.typeBC(i)==6
      
      for i = 1:numel(PARAM.n)
      
          if i==1
              startMatrix = 1;
          else
              startMatrix = 1+sum(nnn(1:i-1));
          end
          endMatrix = sum(nnn(1:i));

          under(2*startMatrix-1:2:2*endMatrix-1) = -8*pi*PARAM.rigidBodyMotion(i);
      
      end
      
  end
  
  %add repulsive force for droplet (or bubble)
  for i = 1:numel(PARAM.panels)
      
      if (PARAM.repulsiveForces(i)==3||PARAM.repulsiveForces(i)==4||PARAM.repulsiveForces(i)==5) && PARAM.blockType(i)==2
      
      if numel(PARAM.panels)>2
          error('Current implementation works only for two block')
      end
      
      %range
      [~,~,~,~,panelRange] = getBlockCoordinates(x,y,PARAM,i);
      startMatrix = getSingRange(panelRange(1),nnn);
      [~,endMatrix] = getSingRange(panelRange(end),nnn);
      
      for k = 1:numel(PARAM.panels)
          
            if k~=i
                
                if PARAM.repulsiveForces(i)==3||PARAM.repulsiveForces(i)==4
      
                    [dfXrep,dfYrep] = disjoiningPressureBlocks(x,y,i,k,PARAM);
                    if PARAM.repulsiveForces(i)==4
                        v(2*startMatrix-1:2:2*endMatrix-1) = v(2*startMatrix-1:2:2*endMatrix-1) + dfXrep;
                    end
                    v(2*startMatrix:2:2*endMatrix) = v(2*startMatrix:2:2*endMatrix) + dfYrep;
                
                elseif PARAM.repulsiveForces(i)==5
                    
                    [dfXrep,dfYrep] = disjoiningPressureBlocksHamaker(x,y,i,k,PARAM);
                    v(2*startMatrix-1:2:2*endMatrix-1) = v(2*startMatrix-1:2:2*endMatrix-1) + dfXrep;
                    v(2*startMatrix:2:2*endMatrix) = v(2*startMatrix:2:2*endMatrix) + dfYrep;
                    
                end
            
            end
            
      end
    
      end
    
  end
  
  %add undelying flow
  if PARAM.addFlow==1
    under(1:2:end-1-addSize) = 8*pi*PARAM.Uunder;
  elseif PARAM.addFlow==2
    
      for i = 1:numel(PARAM.n)
        %range
        [startMatrix,endMatrix] = getSingRange(i,nnn);
        [uExtens,vExtens] = extens_flow(x{i},y{i},PARAM.Ca);
        under(2*startMatrix-1:2:2*endMatrix-1) = under(2*startMatrix-1:2:2*endMatrix-1) + 8*pi*uExtens;
        under(2*startMatrix:2:2*endMatrix) = under(2*startMatrix:2:2*endMatrix) + 8*pi*vExtens;
      end
      
  end
  
  %solve linear system Ay = b, where b = -Uv + under
  b = -U*v - under;

