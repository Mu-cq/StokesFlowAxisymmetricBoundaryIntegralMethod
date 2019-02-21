%integration weight on block

function weight = integrationWeightBlock(x,y,PARAM,nBlock)

error('Look integrationOnLineWeightAxis')

[xHere,yHere,~,~,panelRange,~] = getBlockCoordinates(x,y,PARAM,nBlock);

if PARAM.orderGeometry(panelRange(1))==0 && PARAM.orderVariable(panelRange(1))==0
    
    weight = 2*pi*Ysing(startMatrix:endMatrix)'.*dl;

elseif PARAM.orderGeometry(panelRange(1))==0 && PARAM.orderVariable(panelRange(1))==1
    
    weight = 2*pi*([Ysing(startMatrix:endMatrix-1)'.*dl/2 0] + [0 Ysing(startMatrix+1:endMatrix)'.*dl/2]);
    
elseif PARAM.orderGeometry(panelRange(1))==1 && PARAM.orderVariable(panelRange(1))==1
    
    weight = int_axis_spline_symmetric_weight(xHere,yHere);
    
else
   
    error('Not implemented')
    
end