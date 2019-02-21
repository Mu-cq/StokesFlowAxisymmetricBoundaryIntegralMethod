%replace interface point too close to walls

function [dfX,dfY] = disjoiningPressureBlocksHamaker(x,y,block2,block1,PARAM)

    %panel from which the force is coming
    [xGiveForce,yGiveForce,~,~,panelRange1,dl1] = getBlockCoordinates(x,y,PARAM,block1);
    
    %panel getting the force
    [xReceive,yReceive,nx,ny,panelRange2] = getBlockCoordinates(x,y,PARAM,block2);
    
    if PARAM.orderVariableStokes(panelRange1(1))==0
        
        Xsing1 = (xGiveForce(1:end-1)+xGiveForce(2:end))/2;
        Ysing1 = (yGiveForce(1:end-1)+yGiveForce(2:end))/2;
        
        if PARAM.orderVariableStokes(panelRange2(1))==1
            
            Xsing2 = xReceive;
            Ysing2 = yReceive;
            
        else
            
            error('Not implemented');
            
        end
        
    else
        
        error('Not implemented')
        
    end
    
    %compute distance
    dist = geometryRepulsiveDoubleLayer(Xsing1,Ysing1,Xsing2,Ysing2);
    dist = min(dist')';

    %compute intensity repulsive force
    df = PARAM.coeffRepulsive./dist.^3;
    dfX = df.*nx';
    dfY = df.*ny';
    
%     if PARAM.smoothingRep==1
%     
%         %smooth oscillations
%         dx = diff(xReceive);    dy = diff(yReceive);    dl = sqrt(dx.^2+dy.^2);  l = [0 cumsum(dl)];
%         [~,dfX12] = spaps(l,dfX12,1e-2*max(abs(dfX12)));
%         dfX12 = dfX12';
%         [~,dfY12] = spaps(l,dfY12,1e-2*max(abs(dfY12)));
%         dfY12 = dfY12'.*(dfY12'>0);
%     
%     end

end