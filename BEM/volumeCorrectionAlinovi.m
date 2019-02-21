%deflation on block where A is the influnce matrix and b is known term
%vector

function [Aout,bout] = volumeCorrectionAlinovi(Xsing,Ysing,nnn,x,y,A,b,PARAM)

  if sum(PARAM.blockVolCorr)>0

  %add as many lines as the number of LAGRANGE multiplier
  b = [b; zeros(sum(PARAM.blockVolCorr),1)];
  A = [A; zeros(sum(PARAM.blockVolCorr),2*numel(Xsing))];
  A = [A zeros(2*numel(Xsing)+1,1)];

  countVolCorrBlock = 0;
  for i = 1:numel(PARAM.panels)
      
      %range
      [~,~,~,~,panelRange] = getBlockCoordinates(x,y,PARAM,i);
      startMatrix = getSingRange(panelRange(1),nnn);
      [~,endMatrix] = getSingRange(panelRange(end),nnn);
      
      %index for the first panel of the block
      firstPanel = panelRange(1);
      
      %error('Not working well')
      
      if PARAM.blockVolCorr(i)==1
          
          if PARAM.blockType(i)~=2
          
            error('This is implemented only for a closed interface')
          
          end
          
          countVolCorrBlock = countVolCorrBlock+1;
          
          %integration weights
          weight = integrationOnLineWeightAxis(x{firstPanel},y{firstPanel},PARAM.orderVariable(firstPanel),PARAM.orderGeometry(firstPanel),PARAM.SPlinesType(firstPanel));
          
          %continuity with divergence theorem
          [nx,ny] = computeNormalVector(x{panelRange},y{panelRange},PARAM.orderVariable(panelRange),PARAM.orderGeometry(panelRange),PARAM.SPlinesType(panelRange));
          contDiv = zeros(1,2*numel(weight));
          contDiv(1:2:end-1) = nx.*weight;
          contDiv(2:2:end) = ny.*weight;
          A(2*numel(Xsing)+countVolCorrBlock,2*startMatrix-1:2*endMatrix) = contDiv;
               
          %add column for lagrange multiplier
          A(2*startMatrix-1:2*endMatrix,2*numel(Xsing)+countVolCorrBlock) = contDiv';
          
      end
      
  end
  
  end
  
  Aout = A;
  bout = b;
 
 
 
 
 
 
 
 
 
 
 
 
 