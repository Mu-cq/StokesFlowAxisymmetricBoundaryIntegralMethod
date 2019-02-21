%deflation on block where A is the influnce matrix and b is known term
%vector

function [Aout,bout] = forceFreeStokesAxis(Xsing,Ysing,nnn,x,y,A,b,PARAM)

  countMovingBlock = 0;
  for i = 1:numel(PARAM.panels)
      
      %range
      [xHere,yHere,nxHere,nyHere,panelRange] = getBlockCoordinates(x,y,PARAM,i);
      startMatrix = getSingRange(panelRange(1),nnn);
      [~,endMatrix] = getSingRange(panelRange(end),nnn);
      
      %index for the first panel of the block
      firstPanel = panelRange(1);
      
      if PARAM.blockType(i)==1
          
          countMovingBlock = countMovingBlock+1;
      
          if PARAM.orderGeometry(firstPanel)==0 && PARAM.orderVariable(firstPanel)==0

              %integration in order to impose force free condition
               dl = sqrt(diff(xHere).^2+diff(yHere).^2);
               
               if (PARAM.typeBCstokes(firstPanel)==0||PARAM.typeBCstokes(firstPanel)==1||PARAM.typeBCstokes(firstPanel)==2||PARAM.typeBCstokes(firstPanel)==3||PARAM.typeBCstokes(firstPanel)==4)
                   
                   weight = 2*pi*Ysing(startMatrix:endMatrix)'.*dl;
                   A(2*numel(Xsing)+countMovingBlock,2*startMatrix-1:2:2*endMatrix-1) = weight;
               
               elseif PARAM.typeBCstokes(firstPanel)==5
                   
                   weight = 2*pi*Ysing(startMatrix:endMatrix)'.*dl.*nxHere;
                   A(2*numel(Xsing)+countMovingBlock,2*startMatrix-1:2:2*endMatrix-1) = weight;
                   
                   %tangent stress
                   fT = PARAM.stressBC{firstPanel}.*ones(numel(nxHere),1);
                   weight = -2*pi*Ysing(startMatrix:endMatrix)'.*dl.*nyHere;
                   b(2*numel(Xsing)+countMovingBlock) = b(2*numel(Xsing)+countMovingBlock) + weight*fT;
                   
               else
                   
                   error('Not implemented')
                   
               end
               
               %coefficient of the constant velocity
               A(2*startMatrix+1:2:2*endMatrix-1,2*numel(Xsing)+countMovingBlock) = -8*pi;
               
               %add repulsive forces
               if PARAM.repulsiveForces(i)>0
                   
                  for k = 1:numel(PARAM.blockType)
                      
                      if k~=i
                          
                          [xHere2,yHere2,~,~,panelRange2] = getBlockCoordinates(x,y,PARAM,k);
                          if PARAM.repulsiveForces(i)==1
                              
                              Frep = computeRepulsiveForce2(x,y,i,k,PARAM);
                           
                          elseif PARAM.repulsiveForces(i)==2
                              
                              Frep = computeRepulsiveForce3(xHere,yHere,xHere2,yHere2,PARAM);
                              
                          elseif PARAM.repulsiveForces(i)==3
                              
                              error('Not implemented')
                            
                          elseif PARAM.repulsiveForces(i)==4 || PARAM.repulsiveForces(i)==5
                              
                              if numel(PARAM.panels)>2 || numel(PARAM.blockType(k))==2
                                  
                                  error('Current implementation works only with two block and when second block is a droplet')
                                  
                              end
                              
                              %compute stresses on other block
                              nx = computeNormalVector(x{panelRange2(1)},y{panelRange2(1)},PARAM.orderVariable(panelRange2(1)),PARAM.orderGeometryStokes(panelRange2(1)),PARAM.SPlinesType(panelRange2(1)));
                              [k1,k2] = computeCurvatureSplines(x{panelRange2(1)},y{panelRange2(1)},PARAM.orderVariableStokes(panelRange2(1)));
                              curv = k1+k2;
                              dfX = PARAM.stressBC{panelRange2(1)}*curv.*nx;
                              
                              %compute repulsive forces acting on other block
                              if PARAM.repulsiveForces(i)==4
                                  dfXrep = disjoiningPressureBlocks(x,y,k,i,PARAM);
                              elseif PARAM.repulsiveForces(i)==5
                                  dfXrep = disjoiningPressureBlocksHamaker(x,y,k,i,PARAM);
                              end
                              
                              %integrate
                              weight = integrationOnLineWeightAxis(x{panelRange2(1)},y{panelRange2(1)},PARAM.orderVariableStokes(panelRange2(1)),PARAM.orderGeometryStokes(panelRange2(1)),PARAM.SPlinesType(panelRange2(1)));
                              Frep = -weight*(dfX'+dfXrep);
                              
                          else
                              
                              error('Not implemented')
                           
                          end
                          
                          b(2*numel(Xsing)+countMovingBlock) = b(2*numel(Xsing)+countMovingBlock)+Frep;
                      end
                  end
                   
               end

          elseif PARAM.orderGeometry(panelRange(1))==1 && PARAM.orderVariable(panelRange(1))==0

               error('Not Implemented')

          elseif PARAM.orderGeometry(panelRange(1))==0 && PARAM.orderVariable(panelRange(1))==1
              
              if PARAM.repulsiveForces(i)>0
                  
                  error('Not implemented')
                  
              end
              
              if (PARAM.typeBCstokes(firstPanel)==1||PARAM.typeBCstokes(firstPanel)==2||PARAM.typeBCstokes(firstPanel)==3||PARAM.typeBCstokes(firstPanel)==4)

                   %integration in order to impose force free condition
                   dl = sqrt(diff(xHere).^2+diff(yHere).^2);
                   weight = 2*pi*([Ysing(startMatrix:endMatrix-1)'.*dl/2 0] + [0 Ysing(startMatrix+1:endMatrix)'.*dl/2]);
                   A(2*numel(Xsing)+countMovingBlock,2*startMatrix-1:2:2*endMatrix-1) = weight;
               
               elseif PARAM.typeBCstokes(firstPanel)==5
                   
                   warning('Not validated')
                   weight = 2*pi*([Ysing(startMatrix:endMatrix-1)'.*dl/2 0] + [0 Ysing(startMatrix+1:endMatrix)'.*dl/2]).*nxHere;
                   A(2*numel(Xsing)+countMovingBlock,2*startMatrix-1:2:2*endMatrix-1) = weight;
                   
                   %tangent stress
                   fT = PARAM.stressBC{firstPanel}.*ones(numel(nxHere),1);
                   weight = -2*pi*([Ysing(startMatrix:endMatrix-1)'.*dl/2 0] + [0 Ysing(startMatrix+1:endMatrix)'.*dl/2]).*nyHere;
                   b(2*numel(Xsing)+countMovingBlock) = 2*numel(Xsing)+countMovingBlock + weight*fT;
                   
               else
                   
                  errro('Not implemented')
                   
               end

               %coefficient of the constant velocity
               A(2*startMatrix+1:2:2*endMatrix-1,2*numel(Xsing)+countMovingBlock) = -8*pi;

          elseif PARAM.orderGeometry(panelRange(1))==1 && PARAM.orderVariable(panelRange(1))==1
              
              if PARAM.repulsiveForces(i)>0
                  
                  error('Not implemented')
                  
              end

               error('Not validated')
               weightSPLINE = int_axis_spline_symmetric_weight(xHere,yHere);
               A(2*numel(Xsing)+countMovingBlock,2*startMatrix-1:2:2*endMatrix-1) = weightSPLINE;

          end
      
      end
      
  end
  
  Aout = A;
  bout = b;
 
 
 
 
 
 
 
 
 
 
 
 
 