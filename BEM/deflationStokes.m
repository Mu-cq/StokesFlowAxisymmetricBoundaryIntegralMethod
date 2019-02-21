%deflation on block where A is the influnce matrix and b is known term
%vector

function [Aout,bout] = deflationStokes(Ysing,nnn,x,y,A,b,PARAM)

 %check if deflation is needed
 if sum(PARAM.deflationBlock)>0
      
    %loop over blocks
      for i = 1:numel(PARAM.panels)
          
          if PARAM.deflationBlock(i)==1
      
            % coordinates of the block
            [~,~,nxHere,nyHere,panelRange,dl] = getBlockCoordinates(x,y,PARAM,i);
            
            %index for the first panel of the block
            firstPanel = panelRange(1);
            
            % indeces of the block
            %start and finish
            if i>1
                start = sum(nnn(1:panelRange(1)-1))+1;
            else
                start = 1;
            end
            finish = sum(nnn(1:panelRange(end)));
            
            if PARAM.blockType(i)~=2   % deflation for wall
                
                if PARAM.orderVariable(panelRange(1))==0 && PARAM.orderGeometry(panelRange(1))==0

                    %eigevector as the normal to the element
                    vEigen = zeros(2*numel(dl),1);
                    vEigen(1:2:end-1) = 2*pi*dl.*nxHere.*Ysing(start:finish)';
                    vEigen(2:2:end) = 2*pi*dl.*nyHere.*Ysing(start:finish)';
                    vEigen = vEigen/norm(vEigen);

                    %preconditioning of influnce matrix
                    P = eye(numel(vEigen))-vEigen*vEigen';
                    A(2*start-1:2*finish,:) = P*A(2*start-1:2*finish,:);

                    %preconditioning of rhs
                    b(2*start-1:2*finish) = P*b(2*start-1:2*finish);

                    %add condition in one line, integral of (sigma dot n) = C,
                    %where C is an arbitrary constant
                    if PARAM.typeBCstokes(firstPanel)==5 %normal and tangent coordinates
                        
                        A(2*start-1,:) = 0;
                        A(2*start-1,2*start-1:2:2*finish-1) = 2*pi*Ysing(start:finish)'.*dl;
                        b(2*start-1) = PARAM.deflationConstant(i);
                        
                    elseif (PARAM.typeBCstokes(firstPanel)==0||PARAM.typeBCstokes(firstPanel)==1||PARAM.typeBCstokes(firstPanel)==2||PARAM.typeBCstokes(firstPanel)==3||PARAM.typeBCstokes(firstPanel)==4||PARAM.typeBCstokes(firstPanel)==6) %cartesian coordinates
                        
                        A(2*start-1,:) = 0;
                        A(2*start-1,2*start-1:2:2*finish-1) = 2*pi*nxHere.*Ysing(start:finish)'.*dl;
                        A(2*start-1,2*start:2:2*finish) = 2*pi*nyHere.*Ysing(start:finish)'.*dl;
                        b(2*start-1) = PARAM.deflationConstant(i);
                        
                    else
                        
                        error('Not implemented')
                        
                    end
                    
                elseif PARAM.orderVariable(panelRange(1))==0 && PARAM.orderGeometry(panelRange(1))==1

                    %eigevector as the normal to the element
                    vEigen = zeros(2*numel(dl),1);
                    vEigen(1:2:end-1) = 2*pi*dl.*nxHere.*Ysing(start:finish)';
                    vEigen(2:2:end) = 2*pi*dl.*nyHere.*Ysing(start:finish)';
                    vEigen = vEigen/norm(vEigen);

                    %preconditioning of influnce matrix
                    P = eye(numel(vEigen))-vEigen*vEigen';
                    A(2*start-1:2*finish,:) = P*A(2*start-1:2*finish,:);

                    %preconditioning of rhs
                    b(2*start-1:2*finish) = P*b(2*start-1:2*finish);

                    %add condition in one line, integral of (sigma dot n) = C,
                    %where C is an arbitrary constant
                    if PARAM.typeBCstokes(firstPanel)==5 %normal and tangent coordinates
                        
                        A(2*start-1,:) = 0;
                        A(2*start-1,2*start-1:2:2*finish-1) = 2*pi*Ysing(start:finish)'.*dl;
                        b(2*start-1) = PARAM.deflationConstant(i);
                        
                    elseif (PARAM.typeBCstokes(firstPanel)==1||PARAM.typeBCstokes(firstPanel)==2||PARAM.typeBCstokes(firstPanel)==3||PARAM.typeBCstokes(firstPanel)==4) %cartesian coordinates
                        
                        A(2*start-1,:) = 0;
                        A(2*start-1,2*start-1:2:2*finish-1) = 2*pi*nxHere.*Ysing(start:finish)'.*dl;
                        A(2*start-1,2*start:2:2*finish) = 2*pi*nyHere.*Ysing(start:finish)'.*dl;
                        b(2*start-1) = PARAM.deflationConstant(i);
                        
                    else
                        
                        error('Not implemented')
                        
                    end

                else

                   error('Not implemented')

                end
                
            elseif PARAM.blockType(i)==2   % deflation for bubbles, panel corresponds to block
                
                if PARAM.visc>0.1
                    error('No need for deflation on this block')
                end
                
                if PARAM.orderVariable(panelRange(1))==1 && PARAM.orderGeometry(panelRange(1))==1
                
                    nodes = nnn(panelRange(1));
                    xPanel = x{panelRange(1)};
                    yPanel = y{panelRange(1)};
                    
                    %normal vector
                    nx = nxHere;     ny = nyHere;
                    nn = zeros(1,2*nodes);
                    nn(1:2:end-1) = nx;
                    nn(2:2:end) = ny;
                    nnn = repmat(nn,2*nodes,1);
                    nn = nnn';

                    %integration weights
                    weight = zeros(1,2*nodes);
                    weightSPLINE = int_axis_spline_symmetric_weight(xPanel,yPanel);
                    weight(1:2:end-1) = weightSPLINE;
                    weight(2:2:end) = weightSPLINE;
                    intWeight = repmat(weight,2*nodes,1);

                    %operator for integration
                    INTop = intWeight.*nn.*nnn;
                    A(2*start-1:2*finish,2*start-1:2*finish) = A(2*start-1:2*finish,2*start-1:2*finish) + INTop;

                    %impose volume growth
                    b(2*start-1:2*finish) = b(2*start-1:2*finish) + PARAM.Qsource(i)*diag(nnn);
                
                else
                    
                   error('Not implemented') 
                    
                end
                
            end
        
          end
        
      end
        
 end
 
 Aout = A;
 bout = b;
 
 
 
 
 
 
 
 
 
 
 
 
 