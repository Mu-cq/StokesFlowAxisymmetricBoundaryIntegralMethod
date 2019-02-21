%replace interface point too close to walls

function [dfX12,dfY12,dfX21,dfY21] = disjoiningPressureDoubleLayerDirection2(a1,b1,a2,b2,PARAM)

    %compute distances and direction
    [dist,nx,ny,B1,B2,dl1,dl2] = geometryRepulsiveDoubleLayer(a1,b1,a2,b2);

    %critical distance
    delta = PARAM.repulsiveOn;
    r = dist;

    %compute only when surfaces are close
    activate = r<=delta;

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
    
    %force from 1 to 2
    dfX12 = sum((dfX(:,1:end-1).*dl1.*B1(:,1:end-1)+dfX(:,2:end).*dl1.*B1(:,2:end))/2,2);
    dfY12 = sum((dfY(:,1:end-1).*dl1.*B1(:,1:end-1)+dfY(:,2:end).*dl1.*B1(:,2:end))/2,2);
    %force from 2 to 1
    %dfX21 = -sum(dfX.*dl2'.*B2)';
    %dfY21 = -sum(dfY.*dl2'.*B2)';
    
    %smooth
    dl2 = dl2(1,:);
    l2 = [0 cumsum(dl2)];
    fitobject = fit(l2',dfX12,'poly2');
%     p1 = fitobject.p1;
%     p2 = fitobject.p2;
%     p3 = fitobject.p3;
    
    figure
    plot(fitobject,l2,dfX12)
    grid on
    
%     fFit = pi;
%     figure
%     plot
    

end