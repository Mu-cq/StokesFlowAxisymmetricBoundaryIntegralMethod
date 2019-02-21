%replace interface point too close to walls

function [dfX12,dfY12,dfX21,dfY21] = disjoiningPressureDoubleLayerDirection(a1,b1,a2,b2,PARAM)

    %compute distances and direction
    [dist,nx,ny,B1,B2,dl1,dl2] = geometryRepulsiveDoubleLayer(a1,b1,a2,b2);

    %critical distance
    delta = PARAM.repulsiveOn;
    r = dist;

    %compute only when surfaces are close
    activate = r<=delta;
    
%     if sum(activate)>0
%        here = 1; 
%     end

    %compute intensity repulsive force
    df = PARAM.coeffRepulsive*(exp(delta-r)-1)./(exp(delta)-1).*activate;
    dfX = df.*nx;
    dfY = df.*ny;
    
    %force from 1 to 2
%     dfX12 = sum(dfX,2);
%     dfY12 = sum(dfY,2);
%     %force from 2 to 1
%     dfX21 = -sum(dfX)';
%     dfY21 = -sum(dfY)';
    
    dfX12 = 2*pi*sum(dfX.*[dl1 dl1(:,end)].*B1,2);
    dfY12 = 2*pi*sum(dfY.*[dl1 dl1(:,end)].*B1,2);
    dfX21 = -sum(dfX.*[dl2 dl2(:,end)]'.*B2,2)';
    dfY21 = -sum(dfY.*[dl2 dl2(:,end)]'.*B2,2)';
    %dfX12 = sum((dfX(:,1:end-1).*dl1.*B1(:,1:end-1)+dfX(:,2:end).*dl1.*B1(:,2:end))/2,2);
    %dfY12 = sum((dfY(:,1:end-1).*dl1.*B1(:,1:end-1)+dfY(:,2:end).*dl1.*B1(:,2:end))/2,2);
    %force from 2 to 1
    %dfX21 = -sum(dfX.*dl2'.*B2)';
    %dfY21 = -sum(dfY.*dl2'.*B2)';

end